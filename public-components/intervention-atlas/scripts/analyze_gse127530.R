suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(fgsea)
  library(org.Hs.eg.db)
  library(reactome.db)
  library(jsonlite)
})
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[[1L]])
source(file.path(dirname(script), "bridge_omics_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  index <- match(name, args)
  if (is.na(index) || index == length(args)) default else args[[index + 1L]]
}

contrast_timepoint <- function(contrast_name) {
  labels <- c(three_hr_vs_fast = "3hr", six_hr_vs_fast = "6hr",
              six_hr_vs_three_hr = "6hr_vs_3hr")
  if (length(contrast_name) != 1L || is.na(contrast_name) ||
      !(contrast_name %in% names(labels))) stop("Unknown expression contrast")
  unname(labels[[contrast_name]])
}

data_root <- normalizePath(arg_value("--data-dir", Sys.getenv("NUTRIGENOMICS_DATA_DIR", "data")), mustWork = TRUE)
output_root <- arg_value("--output-dir", file.path(data_root, "..", "outputs", "gse127530-v1"))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
count_dir <- file.path(data_root, "raw", "geo", "GSE127530")
count_file <- file.path(count_dir, "GSE127530_fixed_combinedCounts.txt.gz")
audit_file <- file.path(output_root, "gse127530_counts_audit.json")
if (!file.exists(audit_file)) {
  stop("Run scripts.audit_gse127530_counts first and write its report to the output directory.")
}
count_audit <- jsonlite::fromJSON(audit_file, simplifyVector = FALSE)
if (!(count_audit$status %in% c("identical_payloads", "verified_fixed_cleanup"))) {
  stop("Count file audit did not verify byte-identical payloads.")
}

count_connection <- gzfile(count_file, open = "rt")
counts_table <- read.delim(
  count_connection,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
close(count_connection)
if (ncol(counts_table) != 45L || is.null(rownames(counts_table))) {
  stop("Unexpected GSE127530 table shape or missing gene identifiers.")
}
if (anyDuplicated(rownames(counts_table)) > 0L) {
  stop("Duplicate gene symbols in input; resolve identifiers before analysis.")
}
counts <- as.matrix(counts_table)
storage.mode(counts) <- "numeric"
if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
  stop("Input includes non-finite, negative, or non-integer counts.")
}
if (anyDuplicated(colnames(counts)) > 0L) stop("Duplicate sample labels.")
source_gene_rows <- nrow(counts)
duplicate_gene_symbol_rows <- source_gene_rows - length(unique(rownames(counts)))
if (duplicate_gene_symbol_rows > 0L) {
  counts <- rowsum(counts, group = rownames(counts), reorder = FALSE)
}

sample_pattern <- regexec("^(S[0-9]+)-D([0-9]+)-(Fast|3hr|6hr)_S([0-9]+)$", colnames(counts), ignore.case = TRUE)
parsed <- regmatches(colnames(counts), sample_pattern)
if (any(lengths(parsed) != 5L)) stop("A sample label could not be parsed into person/day/timepoint.")
sample_info <- data.frame(
  sample_label = colnames(counts),
  subject_id = vapply(parsed, function(x) x[[2L]], character(1)),
  study_day = factor(vapply(parsed, function(x) x[[3L]], character(1)), levels = c("1", "2", "3")),
  timepoint = factor(tolower(vapply(parsed, function(x) x[[4L]], character(1))), levels = c("fast", "3hr", "6hr")),
  library_suffix = as.integer(vapply(parsed, function(x) x[[5L]], character(1))),
  stringsAsFactors = FALSE
)
`%||%` <- function(value, fallback) if (is.null(value)) fallback else value
geo_accessions <- vapply(count_audit$sample_design$samples, function(row) row$geo_sample_accession %||% NA_character_, character(1))
if (length(geo_accessions) != ncol(counts) || anyNA(geo_accessions) || anyDuplicated(geo_accessions) > 0L) {
  stop("All count-matrix columns must map one-to-one to GEO GSM sample records before analysis.")
}
sample_info$geo_sample_accession <- geo_accessions
sample_info$subject_id <- factor(sample_info$subject_id)
rownames(sample_info) <- sample_info$sample_label
if (nlevels(sample_info$subject_id) != 5L || any(table(sample_info$subject_id, sample_info$study_day, sample_info$timepoint) != 1L)) {
  stop("Expected one sample per person × study day × timepoint across 5 people, 3 days, and 3 timepoints.")
}
write.table(sample_info, file.path(output_root, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

design <- model.matrix(~ study_day + timepoint, data = sample_info)
if (qr(design)$rank != ncol(design)) stop("The prespecified design matrix is not full rank.")
contrast_matrix <- makeContrasts(
  three_hr_vs_fast = timepoint3hr,
  six_hr_vs_fast = timepoint6hr,
  six_hr_vs_three_hr = timepoint6hr - timepoint3hr,
  levels = design
)

y <- DGEList(counts = counts, genes = data.frame(gene_symbol = rownames(counts)))
keep <- filterByExpr(y, design = design)
if (sum(keep) < 1000L) stop("Too few genes passed edgeR::filterByExpr; inspect input before continuing.")
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
voom_initial <- voom(y, design, plot = FALSE)
block_correlation <- duplicateCorrelation(voom_initial, design, block = sample_info$subject_id)
if (!is.finite(block_correlation$consensus.correlation)) stop("Repeated-measures correlation could not be estimated.")
voom_fit <- voom(
  y,
  design,
  block = sample_info$subject_id,
  correlation = block_correlation$consensus.correlation,
  plot = FALSE
)
fit <- lmFit(
  voom_fit,
  design,
  block = sample_info$subject_id,
  correlation = block_correlation$consensus.correlation
)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit, robust = TRUE)

write_gzip_tsv <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(data, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

gene_tables <- lapply(seq_len(ncol(contrast_matrix)), function(index) {
  contrast_name <- colnames(contrast_matrix)[[index]]
  results <- topTable(fit, coef = index, number = Inf, sort.by = "none", adjust.method = "BH")
  standard_error <- fit$stdev.unscaled[, index] * sqrt(fit$s2.post)
  data.frame(
    feature_type = "gene",
    feature_id = rownames(results),
    feature_label = rownames(results),
    contrast = contrast_name,
    timepoint = contrast_timepoint(contrast_name),
    effect_estimate = results$logFC,
    standard_error = standard_error[rownames(results)],
    confidence_low = results$logFC - qt(.975,fit$df.total)*standard_error[rownames(results)],
    confidence_high = results$logFC + qt(.975,fit$df.total)*standard_error[rownames(results)],
    p_value = results$P.Value,
    q_value = results$adj.P.Val,
    effect_unit = "log2 fold change",
    status = "measured",
    n_people = nlevels(sample_info$subject_id),
    n_samples = nrow(sample_info),
    stringsAsFactors = FALSE
  )
})
gene_results <- do.call(rbind, gene_tables)
gene_results$q_value_contrast <- gene_results$q_value
gene_results$q_value_study_family <- p.adjust(gene_results$p_value,"BH")
gene_results$q_value <- gene_results$q_value_study_family
write_gzip_tsv(gene_results, file.path(output_root, "gene_effects.tsv.gz"))

symbol_to_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(gene_results$feature_id),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "CharacterList"
)
symbol_to_entrez <- vapply(as.list(symbol_to_entrez),function(ids) {
  ids <- unique(as.character(ids));ids <- ids[!is.na(ids)&grepl("^[0-9]+$",ids)]
  if(length(ids)==1L)ids[[1L]] else NA_character_
},character(1))
write.table(data.frame(source_feature=names(symbol_to_entrez),entrez_id=unname(symbol_to_entrez),
  exact_unique_annotation=!is.na(symbol_to_entrez)),file.path(output_root,"feature_crosswalk.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
frozen_reactome <- freeze_reactome(output_root)
reactome_membership <- AnnotationDbi::toTable(reactomeEXTID2PATHID)
gene_columns <- intersect(c("gene_id", "ENTREZID", "entrez_id"), names(reactome_membership))
path_columns <- intersect(c("DB_ID", "path_id", "PATHID", "pathway_id"), names(reactome_membership))
if (!length(gene_columns) || !length(path_columns)) stop("Unable to identify the gene/pathway columns in Reactome.db.")
gene_col <- gene_columns[[1L]]
path_col <- path_columns[[1L]]
reactome_membership <- reactome_membership[grepl("^R-HSA-", reactome_membership[[path_col]]), , drop = FALSE]
pathway_sets <- split(as.character(reactome_membership[[gene_col]]), as.character(reactome_membership[[path_col]]))
pathway_names_table <- AnnotationDbi::toTable(reactomePATHID2NAME)
path_id_columns <- intersect(c("DB_ID", "path_id", "PATHID", "pathway_id"), names(pathway_names_table))
path_name_columns <- intersect(c("path_name", "PATHNAME", "name"), names(pathway_names_table))
if (!length(path_id_columns) || !length(path_name_columns)) stop("Unable to identify Reactome pathway identifier and name columns.")
pathway_names_table <- pathway_names_table[grepl("^R-HSA-", pathway_names_table[[path_id_columns[[1L]]]]), , drop = FALSE]
pathway_names <- setNames(as.character(pathway_names_table[[path_name_columns[[1L]]]]), as.character(pathway_names_table[[path_id_columns[[1L]]]]))

set.seed(20261001L)
pathway_tables <- lapply(seq_len(ncol(contrast_matrix)), function(index) {
  contrast_name <- colnames(contrast_matrix)[[index]]
  statistic <- fit$t[, index]
  entrez <- unname(symbol_to_entrez[rownames(fit$coefficients)])
  valid <- !is.na(entrez) & nzchar(entrez) & is.finite(statistic)
  ranks <- data.frame(entrez = as.character(entrez[valid]), t = statistic[valid], symbol = rownames(fit$coefficients)[valid])
  # Exclude gene aliases with repeated Entrez IDs rather than selecting by observed t.
  ranks <- ranks[!duplicated(ranks$entrez)&!duplicated(ranks$entrez,fromLast=TRUE), ]
  stats <- setNames(ranks$t, ranks$entrez)
  stats <- sort(stats, decreasing = TRUE)
  enrichment <- fgsea(pathways = pathway_sets, stats = stats, minSize = 10L, maxSize = 500L)
  if (!nrow(enrichment)) return(NULL)
  enrichment$pathway_name <- unname(pathway_names[enrichment$pathway])
  enrichment$leading_edge <- vapply(enrichment$leadingEdge, paste, collapse = ";", character(1))
  enrichment$feature_type <- "pathway"
  enrichment$feature_id <- enrichment$pathway
  enrichment$feature_label <- enrichment$pathway_name
  enrichment$contrast <- contrast_name
  enrichment$timepoint <- contrast_timepoint(contrast_name)
  enrichment$effect_estimate <- enrichment$NES
  enrichment$standard_error <- NA_real_
  enrichment$p_value <- enrichment$pval
  enrichment$q_value <- enrichment$padj
  enrichment$effect_unit <- "normalized enrichment score"
  enrichment$status <- "measured"
  enrichment$n_people <- nlevels(sample_info$subject_id)
  enrichment$n_samples <- nrow(sample_info)
  enrichment$metadata_json <- vapply(seq_len(nrow(enrichment)), function(row) {
    as.character(jsonlite::toJSON(list(
      pathway_size = enrichment$size[[row]],
      leading_edge = enrichment$leading_edge[[row]],
      annotation = "Reactome.db from the locked Bioconductor environment",
      priority_theme = grepl("lipid|fatty acid|cholesterol|immune|inflamm|cytokine|circadian", enrichment$pathway_name[[row]], ignore.case = TRUE)
    ), auto_unbox = TRUE))
  }, character(1))
  enrichment[, c("feature_type", "feature_id", "feature_label", "contrast", "timepoint", "effect_estimate", "standard_error", "p_value", "q_value", "effect_unit", "status", "n_people", "n_samples", "metadata_json")]
})
pathway_tables <- pathway_tables[!vapply(pathway_tables, is.null, logical(1))]
if (!length(pathway_tables)) stop("No Reactome pathway enrichment results were generated.")
pathway_results <- do.call(rbind, pathway_tables)
pathway_results$q_value_contrast <- pathway_results$q_value
pathway_results$q_value_study_family <- p.adjust(pathway_results$p_value,"BH")
pathway_results$q_value <- pathway_results$q_value_study_family
write_gzip_tsv(pathway_results, file.path(output_root, "pathway_effects.tsv.gz"))

sample_qc <- data.frame(
  sample_label = colnames(counts),
  library_size_raw = colSums(counts),
  library_size_filtered = y$samples$lib.size,
  tmm_normalization_factor = y$samples$norm.factors,
  subject_id = sample_info$subject_id,
  study_day = sample_info$study_day,
  timepoint = sample_info$timepoint
)
write.table(sample_qc, file.path(output_root, "sample_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

metadata <- list(
  accession = "GSE127530",
  analysis_version = "gse127530-limma-voom-final-v1",
  run_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input_file = basename(count_file),
  input_file_sha256 = count_audit$files[[2L]]$compressed_sha256,
  input_compressed_bytes = count_audit$files[[2L]]$bytes,
  input_gene_rows = source_gene_rows,
  duplicate_gene_symbol_rows_aggregated = duplicate_gene_symbol_rows,
  filtered_gene_rows = sum(keep),
  sample_count = nrow(sample_info),
  independent_people = nlevels(sample_info$subject_id),
  samples_per_person = as.list(table(sample_info$subject_id)),
  timepoints = levels(sample_info$timepoint),
  study_days = levels(sample_info$study_day),
  block_variable = "subject_id",
  duplicate_correlation = unname(block_correlation$consensus.correlation),
  normalization = "edgeR TMM; limma-voom with duplicateCorrelation blocking on person",
  design_formula = "~ study_day + timepoint; person-level repeated-measures correlation blocked by subject_id",
  contrasts = colnames(contrast_matrix),
  multiple_testing = "Primary BH across all genes times three planned contrasts; pathways form a separate family; within-contrast q retained",
  reactome_snapshot = frozen_reactome$annotation,
  pathway_analysis = "fgsea ranked by limma moderated t statistic; Reactome annotation versions recorded through package lock",
  package_versions = as.list(vapply(c("R", "edgeR", "limma", "fgsea", "AnnotationDbi", "org.Hs.eg.db", "reactome.db"), function(package) {
    if (package == "R") as.character(getRversion()) else as.character(packageVersion(package))
  }, character(1))),
  interpretation_limit = paste(
    "Five independent people, with repeated samples across three study days and three timepoints.",
    "The contrast estimates a whole-blood response to a mixed high-fat challenge and cannot identify the causal effect of a single food component.",
    "No individual response or health-risk prediction is produced."
  )
)
jsonlite::write_json(metadata, file.path(output_root, "analysis_metadata.json"), auto_unbox = TRUE, pretty = TRUE)
pdf(file.path(output_root, "voom_mean_variance.pdf"), width = 7, height = 7)
invisible(voom(y, design, block = sample_info$subject_id, correlation = block_correlation$consensus.correlation, plot = TRUE))
dev.off()
cat(jsonlite::toJSON(list(status = "completed", output_root = output_root, genes = sum(keep), people = nlevels(sample_info$subject_id), samples = nrow(sample_info), correlation = unname(block_correlation$consensus.correlation)), auto_unbox = TRUE, pretty = TRUE), "\n")
