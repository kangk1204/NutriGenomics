# Shared, auditable exports for staged intervention analyses.
suppressPackageStartupMessages({ library(limma); library(jsonlite); library(AnnotationDbi); library(org.Hs.eg.db); library(reactome.db); library(fgsea) })
bridge_args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(name, default) { k <- match(name, bridge_args); if(is.na(k) || k == length(bridge_args)) default else bridge_args[[k+1L]] }
write_gzip_tsv <- function(x,p) { con <- gzfile(p,"wt"); on.exit(close(con)); write.table(x,con,sep="\t",quote=FALSE,row.names=FALSE,na="NA") }
file_sha <- function(path) unname(strsplit(system2("sha256sum",shQuote(path),stdout=TRUE)," +")[[1L]][1L])
read_geo_matrix <- function(path) {
  lines <- readLines(gzfile(path),warn=FALSE)
  begin <- match("!series_matrix_table_begin",lines); end <- match("!series_matrix_table_end",lines)
  if(is.na(begin)||is.na(end)) stop("GEO matrix delimiters missing")
  tab <- read.delim(text=paste(lines[(begin+1L):(end-1L)],collapse="\n"),check.names=FALSE,row.names=1,stringsAsFactors=FALSE)
  result <- as.matrix(tab); storage.mode(result)<-"numeric"
  if(any(!is.finite(result))||anyDuplicated(rownames(result))||anyDuplicated(colnames(result))) stop("Invalid matrix")
  result
}
get_characteristic <- function(record, name) {
  x <- unlist(record$characteristics_ch1); pattern <- paste0("^",name,": ")
  selected <- grep(pattern,x,value=TRUE,ignore.case=TRUE)
  if(length(selected)!=1L) stop(paste("Characteristic absent/ambiguous",name,record$sample_accession))
  sub(pattern,"",selected,ignore.case=TRUE)
}
freeze_reactome <- function(output) {
  members <- AnnotationDbi::toTable(reactomeEXTID2PATHID)
  genes <- intersect(c("gene_id","ENTREZID","entrez_id"),names(members))[1]
  paths <- intersect(c("DB_ID","path_id","PATHID","pathway_id"),names(members))[1]
  members <- data.frame(entrez_id=as.character(members[[genes]]),pathway_id=as.character(members[[paths]]),stringsAsFactors=FALSE)
  members <- unique(members[grepl("^R-HSA-",members$pathway_id),]); members <- members[order(members$pathway_id,members$entrez_id),]
  names_table <- AnnotationDbi::toTable(reactomePATHID2NAME)
  id_col <- intersect(c("DB_ID","path_id","PATHID","pathway_id"),names(names_table))[1]
  name_col <- intersect(c("path_name","PATHNAME","name"),names(names_table))[1]
  names_map <- setNames(as.character(names_table[[name_col]]),as.character(names_table[[id_col]]))
  path <- file.path(output,"reactome_membership.tsv")
  write.table(members,path,sep="\t",quote=FALSE,row.names=FALSE)
  annotation <- list(package="reactome.db",package_version=as.character(packageVersion("reactome.db")),
    file_sha256=file_sha(path),database_metadata=as.list(AnnotationDbi::metadata(reactome.db)),
    release_status="Legacy package snapshot frozen as exact gene-set file; not relabeled as current Reactome release")
  jsonlite::write_json(annotation,file.path(output,"reactome_snapshot.json"),auto_unbox=TRUE,pretty=TRUE)
  list(sets=split(members$entrez_id,members$pathway_id),names=names_map,annotation=annotation)
}
export_expression <- function(expression, sample_info, design, contrasts, accession, version, output,
                              entrez_ids, labels, context, input_files, qc_limit, voom_object=NULL) {
  if(qr(design)$rank!=ncol(design)) stop("Design is not full rank")
  fitted <- lmFit(if(is.null(voom_object)) expression else voom_object,design)
  fit <- eBayes(contrasts.fit(fitted,contrasts),robust=TRUE)
  gene_tables <- lapply(seq_len(ncol(contrasts)),function(i) {
    se <- fit$stdev.unscaled[,i]*sqrt(fit$s2.post); critical <- qt(.975,df=fit$df.total)
    data.frame(feature_type="gene",feature_id=rownames(expression),feature_label=labels,
      contrast=colnames(contrasts)[i],tissue=context$tissue,timepoint=context$timepoint[[i]],
      intervention=context$intervention,dose=context$dose[[i]],duration=context$duration[[i]],
      effect_estimate=fit$coefficients[,i],standard_error=se,confidence_low=fit$coefficients[,i]-critical*se,
      confidence_high=fit$coefficients[,i]+critical*se,p_value=fit$p.value[,i],
      q_value=p.adjust(fit$p.value[,i],"BH"),effect_unit="log2 expression difference",status="measured",
      n_people=length(unique(sample_info$subject_id)),n_samples=nrow(sample_info),stringsAsFactors=FALSE)
  })
  gene_results <- do.call(rbind,gene_tables)
  gene_results$q_value_contrast <- gene_results$q_value
  gene_results$q_value_study_family <- p.adjust(gene_results$p_value,"BH")
  gene_results$q_value <- gene_results$q_value_study_family
  write_gzip_tsv(gene_results,file.path(output,"gene_effects.tsv.gz"))
  annotation <- freeze_reactome(output)
  valid <- !is.na(entrez_ids) & nzchar(entrez_ids)
  selected <- which(valid); selected <- selected[!duplicated(entrez_ids[selected])]
  sets <- lapply(annotation$sets,function(ids) which(entrez_ids[selected] %in% ids))
  sets <- sets[lengths(sets)>=10L & lengths(sets)<=500L]
  camera_tables <- list(); fgsea_tables <- list()
  set.seed(20261001)
  for(i in seq_len(ncol(contrasts))) {
    statistics <- fit$t[selected,i]; names(statistics)<-entrez_ids[selected]
    fg <- as.data.frame(fgsea(pathways=annotation$sets,stats=sort(statistics,decreasing=TRUE),minSize=10L,maxSize=500L))
    fg$leading_edge <- vapply(fg$leadingEdge,paste,collapse=";",character(1))
    fgsea_tables[[i]] <- data.frame(feature_type="pathway",feature_id=fg$pathway,feature_label=unname(annotation$names[fg$pathway]),
      contrast=colnames(contrasts)[i],tissue=context$tissue,timepoint=context$timepoint[[i]],intervention=context$intervention,
      dose=context$dose[[i]],duration=context$duration[[i]],effect_estimate=fg$NES,standard_error=NA_real_,p_value=fg$pval,q_value=fg$padj,
      effect_unit="normalized enrichment score",status=ifelse(is.finite(fg$NES)&is.finite(fg$pval),"measured","not_estimable"),n_people=length(unique(sample_info$subject_id)),n_samples=nrow(sample_info),
      metadata_json=vapply(seq_len(nrow(fg)),function(j) as.character(toJSON(list(leading_edge=fg$leading_edge[j],gene_set_sha256=annotation$annotation$file_sha256,method="fgsea"),auto_unbox=TRUE)),character(1)))
    # CAMERA sensitivity uses one contrast value per independent participant.
    # It must not count a participant's several visits/timepoints as independent.
    people<-unique(as.character(sample_info$subject_id))
    nuisance<-which(rowSums(abs(contrasts))==0)
    paired_values<-sapply(people,function(person){
      cells<-which(as.character(sample_info$subject_id)==person)
      coefficients<-contrasts[,i]
      terms<-which(coefficients!=0)
      result<-rep(0,nrow(expression))
      for(term in terms){
        matching<-cells[design[cells,term]==1]
        if(!length(matching))stop("CAMERA person lacks a contrast cell")
        result<-result+coefficients[term]*rowMeans(expression[,matching,drop=FALSE])
      }
      result
    })
    camera_design<-matrix(1,nrow=length(people),ncol=1,dimnames=list(people,"mean_paired_contrast"))
    if(length(nuisance)) {
      paired_nuisance<-t(sapply(people,function(person){
        cells<-which(as.character(sample_info$subject_id)==person); value<-rep(0,length(nuisance))
        for(term in which(contrasts[,i]!=0)) {
          matching<-cells[design[cells,term]==1]
          value<-value+contrasts[term,i]*colMeans(design[matching,nuisance,drop=FALSE])
        }
        value
      }))
      varying<-apply(paired_nuisance,2,sd)>1e-8
      if(any(varying))camera_design<-cbind(camera_design,paired_nuisance[,varying,drop=FALSE])
    }
    if(qr(camera_design)$rank!=ncol(camera_design))stop("CAMERA paired-period design aliased")
    cam <- camera(paired_values[selected,,drop=FALSE],index=sets,design=camera_design,contrast=1,inter.gene.cor=.01,sort=FALSE)
    camera_tables[[i]] <- data.frame(pathway_id=rownames(cam),contrast=colnames(contrasts)[i],cam,row.names=NULL)
  }
  pathways <- do.call(rbind,fgsea_tables); pathways$q_value_contrast<-pathways$q_value
  pathways$q_value_study_family<-p.adjust(pathways$p_value,"BH"); pathways$q_value<-pathways$q_value_study_family
  write_gzip_tsv(pathways,file.path(output,"pathway_effects.tsv.gz"))
  cameras<-do.call(rbind,camera_tables); cameras$FDR_study_family<-p.adjust(cameras$PValue,"BH")
  write_gzip_tsv(cameras,file.path(output,"camera_sensitivity.tsv.gz"))
  sample_info$inclusion_status<-"included"; sample_info$inclusion_reason<-"Explicit metadata mapping; complete person-condition structure validated"
  write.table(sample_info,file.path(output,"sample_metadata.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  quantile_qc<-as.data.frame(t(apply(expression,2,quantile,c(.01,.25,.5,.75,.99))))
  quantile_qc$sample_accession<-colnames(expression)
  write.table(quantile_qc,file.path(output,"sample_qc.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  metadata<-list(accession=accession,analysis_version=version,run_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    script_path=normalizePath(script),script_sha256=file_sha(script),helper_path=normalizePath(file.path(dirname(script),"bridge_omics_helpers.R")),
    helper_sha256=file_sha(file.path(dirname(script),"bridge_omics_helpers.R")),
    independent_people=length(unique(sample_info$subject_id)),sample_count=nrow(sample_info),
    contrasts=colnames(contrasts),design_columns=colnames(design),person_block="Person fixed effects; biological person is the independent unit",
    multiple_testing="BH within contrast retained; primary q_value BH across all gene-by-primary-contrast tests in this study; pathway family separate",
    input_manifest=lapply(input_files,function(p) list(path=normalizePath(p),sha256=file_sha(p),bytes=file.info(p)$size)),
    result_files=c("gene_effects.tsv.gz","pathway_effects.tsv.gz"),
    auxiliary_files=c("sample_metadata.tsv","sample_qc.tsv","camera_sensitivity.tsv.gz","reactome_membership.tsv","reactome_snapshot.json"),
    normalization=context$normalization,qc_limit=qc_limit,interpretation_limit=context$interpretation,
    pathway_analysis="fgsea moderated-t rankings; CAMERA on one paired contrast per independent person at inter.gene.cor=0.01; varying period differences included as covariates with their model degrees of freedom",
    package_versions=as.list(vapply(intersect(c("limma","edgeR","affy","makecdfenv","hgu133plus2cdf","fgsea","org.Hs.eg.db","reactome.db","jsonlite"),loadedNamespaces()),function(p)as.character(packageVersion(p)),character(1))))
  jsonlite::write_json(metadata,file.path(output,"analysis_metadata.json"),auto_unbox=TRUE,pretty=TRUE)
  saveRDS(list(expression=expression,sample_info=sample_info,design=design,contrasts=contrasts),file.path(output,"reanalysis_inputs.rds"),compress=FALSE)
  cat(toJSON(list(status="completed",accession=accession,genes=nrow(expression),people=metadata$independent_people,
    samples=nrow(sample_info),contrasts=ncol(contrasts),output=output),auto_unbox=TRUE),"\n")
}
