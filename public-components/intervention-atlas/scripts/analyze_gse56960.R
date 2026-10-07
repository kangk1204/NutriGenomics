script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument)) {
  project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_argument[[1L]])), ".."))
  r_library <- Sys.getenv(
    "NUTRIGENOMICS_R_LIB",
    file.path(project_root, "..", "toolchain", "r-library")
  )
  if (dir.exists(r_library)) .libPaths(c(normalizePath(r_library), .libPaths()))
}
Sys.setenv(R_THREADS = "1")

suppressPackageStartupMessages({
  library(affy)
  library(limma)
  library(AnnotationDbi)
  library(hgu219.db)
  library(fgsea)
  library(reactome.db)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  index <- match(name, args)
  if (is.na(index) || index == length(args)) default else args[[index + 1L]]
}

data_root <- normalizePath(arg_value("--data-dir", Sys.getenv("NUTRIGENOMICS_DATA_DIR", "data")), mustWork = TRUE)
output_root <- arg_value("--output-dir", file.path(data_root, "..", "outputs", "gse56960-v1"))
sample_sheet_path <- arg_value(
  "--sample-sheet",
  file.path(data_root, "..", "outputs", "geo-design-audit", "GSE56960", "sample_design.tsv")
)
cel_root <- file.path(data_root, "processed", "geo", "GSE56960", "raw_files")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(sample_sheet_path)) stop("Run scripts.audit_geo_design before analyzing GSE56960.")
sample_info <- read.delim(sample_sheet_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("sample_accession", "title", "subject_candidate", "characteristics_json")
if (!all(required_columns %in% names(sample_info))) stop("GSE56960 sample sheet is missing required columns.")
if (nrow(sample_info) != 168L || anyDuplicated(sample_info$sample_accession)) {
  stop("Expected 168 unique GEO samples for GSE56960.")
}

title_pattern <- paste0(
  "^Subject\\s+([^,]+),\\s*Meal\\s+([A-Z])\\s*\\(([0-9]+)\\s*kcal\\),\\s*",
  "(Fasting|Postprandial)\\s*\\(([0-9]+)h\\)$"
)
parsed_title <- regmatches(sample_info$title, regexec(title_pattern, sample_info$title, ignore.case = TRUE))
if (any(lengths(parsed_title) != 6L)) stop("Could not parse all titles into subject, meal, dose, and timepoint.")
parsed <- do.call(rbind, lapply(parsed_title, function(value) value[-1L]))
sample_info$subject_id <- parsed[, 1L]
if (!identical(as.character(sample_info$subject_candidate), sample_info$subject_id)) {
  stop("Parsed subject identifiers disagree with the metadata audit candidates.")
}
sample_info$meal <- toupper(parsed[, 2L])
sample_info$dose_kcal <- as.integer(parsed[, 3L])
sample_info$timepoint <- paste0(as.integer(parsed[, 5L]), "h")
sample_info$status_group <- ifelse(
  grepl("metabolic status: normal weight", sample_info$characteristics_json, ignore.case = TRUE),
  "NormalWeight",
  ifelse(grepl("metabolic status: obese", sample_info$characteristics_json, ignore.case = TRUE), "Obese", NA_character_)
)
if (anyNA(sample_info$status_group)) stop("Metabolic status is missing or unrecognized for at least one sample.")
if (!identical(sort(unique(sample_info$dose_kcal)), c(500L, 1000L, 1500L))) stop("Unexpected meal doses.")
if (!identical(sort(unique(sample_info$timepoint)), c("0h", "2h", "4h", "6h"))) stop("Unexpected timepoints.")
if (length(unique(sample_info$subject_id)) != 14L) stop("Expected 14 independent participants.")
status_by_subject <- tapply(sample_info$status_group, sample_info$subject_id, function(x) length(unique(x)))
if (any(status_by_subject != 1L)) stop("Metabolic status is inconsistent within a subject.")
if (!all(table(factor(sample_info$status_group, levels = c("NormalWeight", "Obese"))) == 84L)) {
  stop("Expected seven participants and 84 arrays per status group.")
}
cell_key <- paste(sample_info$subject_id, sample_info$meal, sample_info$timepoint, sep = "|")
if (length(unique(cell_key)) != nrow(sample_info)) stop("Duplicate subject × meal × timepoint records were found.")
if (any(table(sample_info$subject_id, sample_info$meal, sample_info$timepoint) != 1L)) {
  stop("Expected one sample per subject × meal × timepoint cell.")
}

cel_gz <- list.files(cel_root, pattern = "(?i)\\.cel\\.gz$", recursive = TRUE, full.names = TRUE)
if (!length(cel_gz)) stop("No compressed CEL files found in the extracted GSE56960 archive.")
gsm_from_file <- sub("_.*$", "", basename(cel_gz))
if (anyDuplicated(gsm_from_file)) stop("Multiple CEL files map to one GSM accession.")
file_index <- match(sample_info$sample_accession, gsm_from_file)
if (anyNA(file_index) || length(cel_gz) != nrow(sample_info)) {
  stop("CEL files and GEO samples must form a complete one-to-one mapping.")
}
cel_gz <- cel_gz[file_index]

inflate_gzip <- function(source_path, destination_path) {
  input <- gzfile(source_path, open = "rb")
  output <- file(destination_path, open = "wb")
  on.exit(close(input), add = TRUE)
  on.exit(close(output), add = TRUE)
  repeat {
    block <- readBin(input, what = "raw", n = 1024L * 1024L)
    if (!length(block)) break
    writeBin(block, output)
  }
}

scratch <- tempfile("gse56960-cel-")
dir.create(scratch, recursive = TRUE)
cel_paths <- file.path(scratch, paste0(sample_info$sample_accession, ".CEL"))
for (index in seq_along(cel_gz)) inflate_gzip(cel_gz[[index]], cel_paths[[index]])

batch <- ReadAffy(filenames = cel_paths, sampleNames = sample_info$sample_accession, verbose = FALSE)
expression_set <- rma(batch, verbose = TRUE)
unlink(scratch, recursive = TRUE)
rm(batch)
invisible(gc())
probe_expression <- exprs(expression_set)
if (!identical(colnames(probe_expression), sample_info$sample_accession)) {
  stop("RMA columns do not match the validated GEO sample order.")
}
probe_ids <- rownames(probe_expression)
probe_symbols <- AnnotationDbi::mapIds(
  hgu219.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "CharacterList"
)
probe_entrez <- AnnotationDbi::mapIds(
  hgu219.db, keys = probe_ids, column = "ENTREZID", keytype = "PROBEID", multiVals = "CharacterList"
)
unique_annotation <- function(values) {
  values <- unique(as.character(values));values <- values[!is.na(values)&nzchar(values)]
  if(length(values)==1L)values[[1L]] else NA_character_
}
probe_symbols <- vapply(as.list(probe_symbols),unique_annotation,character(1))
probe_entrez <- vapply(as.list(probe_entrez),unique_annotation,character(1))
valid_probe <- !is.na(probe_symbols)&!is.na(probe_entrez)&grepl("^[0-9]+$",probe_entrez)
write.table(data.frame(source_feature=probe_ids,source_symbol=unname(probe_symbols),entrez_id=unname(probe_entrez),included=valid_probe,
  mapping_basis="One unambiguous symbol and one Entrez identifier; ambiguous annotations excluded"),file.path(output_root,"feature_crosswalk.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
if(sum(valid_probe)<10000L)stop("Too few uniquely annotated U219 probes")
gene_expression <- avereps(probe_expression[valid_probe,,drop=FALSE],ID=probe_entrez[valid_probe])
gene_entrez <- setNames(rownames(gene_expression),rownames(gene_expression))
gene_labels <- tapply(probe_symbols[valid_probe],probe_entrez[valid_probe],unique_annotation)

sample_info$status_group <- factor(sample_info$status_group, levels = c("NormalWeight", "Obese"))
sample_info$dose_kcal <- factor(sample_info$dose_kcal, levels = c(500, 1000, 1500))
sample_info$timepoint <- factor(sample_info$timepoint, levels = c("0h", "2h", "4h", "6h"))
sample_info$subject_id <- factor(sample_info$subject_id)
sample_info$condition_id <- paste(
  sample_info$status_group, paste0(sample_info$dose_kcal, "kcal"), sample_info$timepoint, sep = "_"
)
condition_levels <- sort(unique(sample_info$condition_id))
sample_info$condition_id <- factor(sample_info$condition_id, levels = condition_levels)
design <- model.matrix(~ 0 + condition_id, data = sample_info)
colnames(design) <- condition_levels
if (qr(design)$rank != ncol(design)) stop("The prespecified condition design matrix is not full rank.")

contrast_expressions <- character()
contrast_context <- list()
for (status_group in levels(sample_info$status_group)) {
  for (dose_kcal in levels(sample_info$dose_kcal)) {
    baseline <- paste(status_group, paste0(dose_kcal, "kcal"), "0h", sep = "_")
    for (post_time in c("2h", "4h", "6h")) {
      post <- paste(status_group, paste0(dose_kcal, "kcal"), post_time, sep = "_")
      contrast_name <- paste("post", status_group, dose_kcal, post_time, "vs_fast", sep = "_")
      contrast_expressions[[contrast_name]] <- paste(post, "-", baseline)
      contrast_context[[contrast_name]] <- list(
        status_group = status_group,
        dose_kcal = as.integer(dose_kcal),
        timepoint = post_time,
        n_people = length(unique(sample_info$subject_id[sample_info$status_group == status_group])),
        n_samples = 2L * length(unique(sample_info$subject_id[sample_info$status_group == status_group]))
      )
    }
  }
}
if (arg_value("--bridge-contrasts", "false") == "true") {
  contrast_expressions <- character(); contrast_context <- list()
  for (status_group in levels(sample_info$status_group)) for(dose_kcal in c("1000","1500")) for(post_time in c("2h","4h","6h")) {
    prefix <- paste0(status_group,"_")
    contrast_name <- paste("dose",status_group,dose_kcal,"vs500",post_time,"change",sep="_")
    contrast_expressions[[contrast_name]] <- paste0("(",prefix,dose_kcal,"kcal_",post_time," - ",prefix,dose_kcal,"kcal_0h) - (",prefix,"500kcal_",post_time," - ",prefix,"500kcal_0h)")
    contrast_context[[contrast_name]] <- list(status_group=status_group,dose_kcal=as.integer(dose_kcal),timepoint=post_time,
      n_people=length(unique(sample_info$subject_id[sample_info$status_group==status_group])),
      n_samples=4L*length(unique(sample_info$subject_id[sample_info$status_group==status_group])))
  }
}
contrast_matrix <- makeContrasts(contrasts = contrast_expressions, levels = design)

correlation <- duplicateCorrelation(gene_expression, design, block = sample_info$subject_id)
if (!is.finite(correlation$consensus.correlation)) stop("Repeated-measures correlation could not be estimated.")
fit <- lmFit(
  gene_expression, design, block = sample_info$subject_id, correlation = correlation$consensus.correlation
)
fit <- eBayes(contrasts.fit(fit, contrast_matrix), robust = TRUE)

write_gzip_tsv <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(data, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

gene_tables <- lapply(seq_len(ncol(contrast_matrix)), function(index) {
  contrast_name <- colnames(contrast_matrix)[[index]]
  context <- contrast_context[[contrast_name]]
  results <- topTable(fit, coef = index, number = Inf, sort.by = "none", adjust.method = "BH")
  standard_error <- fit$stdev.unscaled[, index] * sqrt(fit$s2.post)
  data.frame(
    feature_type = "gene",
    feature_id = rownames(results),
    feature_label = unname(gene_labels[rownames(results)]),
    contrast = contrast_name,
    tissue = "whole blood cells",
    timepoint = context$timepoint,
    intervention = "high-fat test meal challenge",
    dose = paste0(context$dose_kcal, if(arg_value("--bridge-contrasts","false")=="true") " versus 500 kcal" else " kcal"),
    duration = "acute",
    effect_estimate = results$logFC,
    standard_error = standard_error[rownames(results)],
    confidence_low = results$logFC - qt(.975,fit$df.total)*standard_error[rownames(results)],
    confidence_high = results$logFC + qt(.975,fit$df.total)*standard_error[rownames(results)],
    p_value = results$P.Value,
    q_value = results$adj.P.Val,
    effect_unit = "log2 expression difference",
    status = "measured",
    n_people = context$n_people,
    n_samples = nrow(sample_info),
    metadata_json = vapply(seq_len(nrow(results)), function(row) {
      as.character(toJSON(list(status_group = context$status_group, dose_kcal = context$dose_kcal,
        contrast_assay_count=context$n_samples,model_assay_count=nrow(sample_info)), auto_unbox = TRUE))
    }, character(1)),
    stringsAsFactors = FALSE
  )
})
gene_results <- do.call(rbind, gene_tables)
if(arg_value("--bridge-contrasts","false")=="true") {
  gene_results$q_value_contrast<-gene_results$q_value
  gene_results$q_value_study_family<-p.adjust(gene_results$p_value,"BH")
  gene_results$q_value<-gene_results$q_value_study_family
}
write_gzip_tsv(gene_results, file.path(output_root, "gene_effects.tsv.gz"))

reactome_membership <- AnnotationDbi::toTable(reactomeEXTID2PATHID)
gene_columns <- intersect(c("gene_id", "ENTREZID", "entrez_id"), names(reactome_membership))
path_columns <- intersect(c("DB_ID", "path_id", "PATHID", "pathway_id"), names(reactome_membership))
if (!length(gene_columns) || !length(path_columns)) stop("Unable to identify Reactome gene/pathway columns.")
gene_column <- gene_columns[[1L]]
path_column <- path_columns[[1L]]
reactome_membership <- reactome_membership[grepl("^R-HSA-", reactome_membership[[path_column]]), , drop = FALSE]
pathway_sets <- split(as.character(reactome_membership[[gene_column]]), as.character(reactome_membership[[path_column]]))
pathway_names_table <- AnnotationDbi::toTable(reactomePATHID2NAME)
path_id_columns <- intersect(c("DB_ID", "path_id", "PATHID", "pathway_id"), names(pathway_names_table))
path_name_columns <- intersect(c("path_name", "PATHNAME", "name"), names(pathway_names_table))
if (!length(path_id_columns) || !length(path_name_columns)) stop("Unable to identify Reactome pathway names.")
pathway_names_table <- pathway_names_table[grepl("^R-HSA-", pathway_names_table[[path_id_columns[[1L]]]]), , drop = FALSE]
pathway_names <- setNames(
  as.character(pathway_names_table[[path_name_columns[[1L]]]]),
  as.character(pathway_names_table[[path_id_columns[[1L]]]])
)

set.seed(20261001L)
pathway_tables <- lapply(seq_len(ncol(contrast_matrix)), function(index) {
  contrast_name <- colnames(contrast_matrix)[[index]]
  context <- contrast_context[[contrast_name]]
  entrez <- unname(gene_entrez[rownames(fit$coefficients)])
  statistic <- fit$t[, index]
  valid <- !is.na(entrez) & nzchar(entrez) & is.finite(statistic)
  ranks <- data.frame(entrez = as.character(entrez[valid]), t = statistic[valid], symbol = rownames(fit$coefficients)[valid])
  ranks <- ranks[order(-abs(ranks$t), ranks$entrez, ranks$symbol), ]
  ranks <- ranks[!duplicated(ranks$entrez), ]
  stats <- sort(setNames(ranks$t, ranks$entrez), decreasing = TRUE)
  enrichment <- fgsea(pathways = pathway_sets, stats = stats, minSize = 10L, maxSize = 500L)
  if (!nrow(enrichment)) return(NULL)
  enrichment$pathway_name <- unname(pathway_names[enrichment$pathway])
  enrichment$leading_edge <- vapply(enrichment$leadingEdge, paste, collapse = ";", character(1))
  data.frame(
    feature_type = "pathway",
    feature_id = enrichment$pathway,
    feature_label = enrichment$pathway_name,
    contrast = contrast_name,
    tissue = "whole blood cells",
    timepoint = context$timepoint,
    intervention = "high-fat test meal challenge",
    dose = paste0(context$dose_kcal, if(arg_value("--bridge-contrasts","false")=="true") " versus 500 kcal" else " kcal"),
    duration = "acute",
    effect_estimate = enrichment$NES,
    standard_error = NA_real_,
    p_value = enrichment$pval,
    q_value = enrichment$padj,
    effect_unit = "normalized enrichment score",
    status = ifelse(is.finite(enrichment$NES)&is.finite(enrichment$pval),"measured","not_estimable"),
    n_people = context$n_people,
    n_samples = nrow(sample_info),
    metadata_json = vapply(seq_len(nrow(enrichment)), function(row) {
      as.character(toJSON(list(
        status_group = context$status_group,
        pathway_size = enrichment$size[[row]],
        leading_edge = enrichment$leading_edge[[row]],
        annotation = "Reactome.db from the locked Bioconductor environment"
      ), auto_unbox = TRUE))
    }, character(1)),
    stringsAsFactors = FALSE
  )
})
pathway_tables <- pathway_tables[!vapply(pathway_tables, is.null, logical(1))]
if (!length(pathway_tables)) stop("No Reactome pathway results were generated.")
pathway_results <- do.call(rbind, pathway_tables)
if(arg_value("--bridge-contrasts","false")=="true") {
  pathway_results$q_value_contrast<-pathway_results$q_value
  pathway_results$q_value_study_family<-p.adjust(pathway_results$p_value,"BH")
  pathway_results$q_value<-pathway_results$q_value_study_family
  members_file<-file.path(output_root,"reactome_membership.tsv")
  write.table(reactome_membership,members_file,sep="\t",quote=FALSE,row.names=FALSE)
  snapshot<-list(package="reactome.db",package_version=as.character(packageVersion("reactome.db")),
    gene_set_sha256=strsplit(system2("sha256sum",shQuote(members_file),stdout=TRUE)," +")[[1]][1],
    release_status="Frozen legacy package snapshot; not relabeled as a current Reactome release")
  write_json(snapshot,file.path(output_root,"reactome_snapshot.json"),auto_unbox=TRUE,pretty=TRUE)
  entrez_ids<-unname(gene_entrez[rownames(gene_expression)])
  selected<-which(!is.na(entrez_ids)&nzchar(entrez_ids));selected<-selected[!duplicated(entrez_ids[selected])]
  indices<-lapply(pathway_sets,function(ids)which(entrez_ids[selected]%in%ids));indices<-indices[lengths(indices)>=10&lengths(indices)<=500]
  camera_tables<-lapply(seq_len(ncol(contrast_matrix)),function(i){
    context<-contrast_context[[colnames(contrast_matrix)[i]]]
    participants<-unique(as.character(sample_info$subject_id[sample_info$status_group==context$status_group]))
    delta_expression<-sapply(participants,function(person){
      value<-rep(0,nrow(gene_expression))
      for(term in which(contrast_matrix[,i]!=0)) {
        cell<-which(as.character(sample_info$subject_id)==person & design[,term]==1)
        if(length(cell)!=1L)stop("Independent CAMERA contrast cell absent")
        value<-value+contrast_matrix[term,i]*gene_expression[,cell]
      }
      value
    })
    camera_fit<-camera(delta_expression[selected,,drop=FALSE],index=indices,design=matrix(1,length(participants),1),contrast=1,inter.gene.cor=.01,sort=FALSE)
    data.frame(pathway_id=rownames(camera_fit),contrast=colnames(contrast_matrix)[i],camera_fit,row.names=NULL)
  })
  camera_results<-do.call(rbind,camera_tables);camera_results$FDR_study_family<-p.adjust(camera_results$PValue,"BH")
  write_gzip_tsv(camera_results,file.path(output_root,"camera_sensitivity.tsv.gz"))
}
write_gzip_tsv(pathway_results, file.path(output_root, "pathway_effects.tsv.gz"))

write.table(sample_info, file.path(output_root, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sample_qc <- data.frame(
  sample_accession = colnames(gene_expression),
  median_log2_expression = apply(gene_expression, 2L, median, na.rm = TRUE),
  mean_log2_expression = colMeans(gene_expression, na.rm = TRUE),
  sd_log2_expression = apply(gene_expression, 2L, sd, na.rm = TRUE),
  subject_id = sample_info$subject_id,
  status_group = sample_info$status_group,
  dose_kcal = sample_info$dose_kcal,
  timepoint = sample_info$timepoint
)
write.table(sample_qc, file.path(output_root, "sample_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pdf(file.path(output_root, "rma_qc.pdf"), width = 10, height = 7)
boxplot(as.data.frame(gene_expression), outline = FALSE, las = 2, cex.axis = 0.35, main = "GSE56960 RMA expression")
plot(density(gene_expression[, 1L]), main = "RMA expression distributions", xlab = "log2 intensity")
for (sample_index in 2:ncol(gene_expression)) lines(density(gene_expression[, sample_index]), col = rgb(0, 0, 0, 0.12))
dev.off()

`%||%` <- function(value, fallback) if (is.null(value)) fallback else value
sha256_file <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  if (!length(output) || (attr(output, "status") %||% 0L) != 0L) stop(paste("SHA-256 failed for", path))
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
script_path <- normalizePath(file.path("scripts", "analyze_gse56960.R"), mustWork = TRUE)
lock_path <- normalizePath(file.path(data_root, "..", "toolchain", "r-environment-linux-64.lock.txt"), mustWork = TRUE)
overlay_lock_path <- normalizePath(file.path(data_root, "..", "toolchain", "r-overlay-lock.txt"), mustWork = TRUE)
input_manifest <- lapply(
  c(
    cel_gz,
    sample_sheet_path,
    file.path(data_root, "raw", "geo", "download_manifest.json"),
    file.path(data_root, "processed", "geo", "GSE56960", "raw_files", "extraction_manifest.json"),
    lock_path,
    overlay_lock_path
  ),
  function(path) {
  list(path = normalizePath(path, mustWork = TRUE), bytes = file.info(path)$size, sha256 = sha256_file(path))
  }
)
metadata <- list(
  accession = "GSE56960",
  analysis_version = arg_value("--analysis-version",if(arg_value("--bridge-contrasts","false")=="true") "gse56960-rma-limma-dose-change-v2" else "gse56960-rma-limma-v1"),
  run_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script_path = script_path,
  script_sha256 = sha256_file(script_path),
  environment_lock_path = lock_path,
  environment_lock_sha256 = sha256_file(lock_path),
  overlay_lock_path = overlay_lock_path,
  overlay_lock_sha256 = sha256_file(overlay_lock_path),
  input_manifest = input_manifest,
  sample_count = nrow(sample_info),
  independent_people = length(unique(sample_info$subject_id)),
  people_per_status = as.list(table(sample_info$status_group) / 12L),
  samples_per_person = as.list(table(sample_info$subject_id)),
  meal_doses_kcal = as.list(levels(sample_info$dose_kcal)),
  timepoints = as.list(levels(sample_info$timepoint)),
  repeated_measures_correlation = correlation$consensus.correlation,
  normalization = "affy::rma background correction, quantile normalization, and probe summarization",
  preprocessing_threads = Sys.getenv("R_THREADS"),
  gene_level_aggregation = "Unambiguous HGU219 probe-to-Entrez mapping; ambiguous targets excluded; limma::avereps within Entrez gene",
  design = "condition means by metabolic status × meal kcal × timepoint; subject-blocked duplicateCorrelation",
  contrast_count = ncol(contrast_matrix),
  fgsea_random_seed = 20261001L,
  multiple_testing = if(arg_value("--bridge-contrasts","false")=="true") "Primary BH across all genes times 12 dose-change contrasts; pathway family separate; within-contrast q retained" else "BH within contrast",
  result_files = list("gene_effects.tsv.gz", "pathway_effects.tsv.gz"),
  auxiliary_files = list("sample_metadata.tsv", "sample_qc.tsv", "rma_qc.pdf"),
  effect_scope = "cohort-level response to a mixed high-fat challenge; not component-specific and not an individual prediction",
  software_versions = list(
    R = as.character(getRversion()),
    affy = as.character(packageVersion("affy")),
    preprocessCore = as.character(packageVersion("preprocessCore")),
    hgu219.db = as.character(packageVersion("hgu219.db")),
    limma = as.character(packageVersion("limma")),
    reactome.db = as.character(packageVersion("reactome.db"))
  )
)
write_json(metadata, file.path(output_root, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(list(
  status = "completed",
  accession = "GSE56960",
  n_samples = nrow(sample_info),
  n_people = length(unique(sample_info$subject_id)),
  genes = nrow(gene_expression),
  gene_effects = nrow(gene_results),
  pathway_effects = nrow(pathway_results),
  contrasts = ncol(contrast_matrix),
  repeated_correlation = correlation$consensus.correlation
), auto_unbox = TRUE, pretty = TRUE), "\n")
