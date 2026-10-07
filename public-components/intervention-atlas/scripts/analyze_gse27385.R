script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]); source(file.path(dirname(script),"bridge_omics_helpers.R"))
root<-arg_value("--root",NULL)
if(is.null(root)||!nzchar(root)) stop("Explicit --root is required")
raw<-file.path(root,"raw","GSE27385"); output<-arg_value("--output-dir",file.path(root,"outputs","gse27385-processed-v1"))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
records<-fromJSON(file.path(raw,"sample_metadata.json"),simplifyVector=FALSE)
titles<-vapply(records,function(x)x$title[[1]],character(1))
parsed<-regmatches(titles,regexec("^HPBMC_([0-9]+)([ABC])_(FISHOIL|FIBRATE|PLACEBO)$",titles))
if(length(records)!=33L || any(lengths(parsed)!=4L)) stop("Unexpected sample title")
sample_info<-data.frame(sample_accession=vapply(records,function(x)x$sample_accession,character(1)),
 subject_id=vapply(parsed,function(x)x[2],character(1)),period=vapply(parsed,function(x)x[3],character(1)),
 arm=vapply(parsed,function(x)x[4],character(1)))
if(length(unique(sample_info$subject_id))!=11L || any(table(sample_info$subject_id,sample_info$arm)!=1L)) stop("Incomplete person-arm mapping")
expression<-read_geo_matrix(file.path(raw,"GSE27385_series_matrix.txt.gz"))
expression<-expression[,sample_info$sample_accession,drop=FALSE]
raw_mode<-arg_value("--raw-cel","false")=="true"
if(raw_mode) {
 suppressPackageStartupMessages({library(affy);library(hgu133plus2cdf)})
 cel_root<-file.path(root,"processed","GSE27385","raw_files");dir.create(cel_root,recursive=TRUE,showWarnings=FALSE)
 archive<-file.path(raw,"GSE27385_RAW.tar");members<-untar(archive,list=TRUE)
 if(any(grepl("(^/|(^|/)\\.\\.(/|$))",members)))stop("Unsafe archive paths")
 untar(archive,exdir=cel_root)
 compressed<-list.files(cel_root,pattern="[.]CEL[.]gz$",full.names=TRUE,recursive=TRUE)
 gsm<-sub("_.*$","",basename(compressed));matching<-match(sample_info$sample_accession,gsm)
 if(length(compressed)!=33L||anyNA(matching)||anyDuplicated(gsm))stop("Raw CEL and GSM mapping not one-to-one")
 cel_files<-sub("[.]gz$","",compressed[matching])
 for(i in seq_along(cel_files))if(!file.exists(cel_files[i])) {
  input<-gzfile(compressed[matching][i],"rb");destination<-file(cel_files[i],"wb")
  repeat {block<-readBin(input,"raw",1024L*1024L);if(!length(block))break;writeBin(block,destination)}
  close(input);close(destination)
 }
 batch<-ReadAffy(filenames=cel_files,sampleNames=sample_info$sample_accession)
 write.table(data.frame(sample_accession=sample_info$sample_accession,mean_pm=colMeans(pm(batch)),median_pm=apply(pm(batch),2,median)),
  file.path(output,"raw_cel_qc.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
 expression<-exprs(rma(batch,verbose=TRUE));rm(batch);gc()
}
# Read the exact depositor GPL570 annotation embedded in the family SOFT.
soft<-readLines(gzfile(file.path(raw,"GSE27385_family.soft.gz")),warn=FALSE)
start<-match("!platform_table_begin",soft); end<-match("!platform_table_end",soft)
if(is.na(start)||is.na(end)) stop("Original platform annotation missing")
platform<-read.delim(text=paste(soft[(start+1):(end-1)],collapse="\n"),check.names=FALSE,stringsAsFactors=FALSE,quote="",comment.char="")
entrez_col<-intersect(c("ENTREZ_GENE_ID","Entrez Gene","ENTREZ_GENE"),names(platform))[1]
if(is.na(entrez_col)) stop("Entrez annotation not found")
source_entrez<-as.character(platform[[entrez_col]][match(rownames(expression),platform$ID)])
valid<-!is.na(source_entrez)&grepl("^[0-9]+$",source_entrez)
write.table(data.frame(source_id=rownames(expression),source_entrez=source_entrez,included=valid,
 reason=ifelse(valid,"one Entrez identifier","Missing or multiple target genes; excluded from gene-level inference")),
 file.path(output,"feature_crosswalk.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
expression<-avereps(expression[valid,,drop=FALSE],ID=source_entrez[valid]); entrez<-rownames(expression)
labels<-unname(mapIds(org.Hs.eg.db,keys=entrez,column="SYMBOL",keytype="ENTREZID",multiVals="first")); labels[is.na(labels)]<-entrez[is.na(labels)]
sample_info$subject_id<-factor(sample_info$subject_id); sample_info$arm<-factor(sample_info$arm); sample_info$period<-factor(sample_info$period)
design<-model.matrix(~0+arm+subject_id+period,sample_info);colnames(design)<-sub("^arm","",colnames(design))
contrasts<-makeContrasts(fishoil_vs_placebo_6weeks=FISHOIL-PLACEBO,fenofibrate_vs_placebo_6weeks=FIBRATE-PLACEBO,levels=design)
export_expression(expression,sample_info,design,contrasts,"GSE27385",if(raw_mode)"gse27385-rma-limma-paired-period-v2" else "gse27385-limma-paired-period-processed-v1",output,entrez,labels,
 context=list(tissue="PBMC",timepoint=as.list(rep("end of 6-week period",2)),dose=list("3.7 g/day n-3 LCPUFA (EPA 1.7; DHA 1.2 g/day)","200 mg/day fenofibrate"),duration=as.list(rep("6 weeks",2)),
 intervention="Fish oil or fenofibrate versus placebo; randomized crossover",normalization=if(raw_mode)"Independent raw CEL RMA; unambiguous GPL570 Entrez probes averaged per gene; person and period fixed effects" else "Depositor RMA matrix; unambiguous GPL570 Entrez probes averaged per gene; person and period fixed effects",
 interpretation="11 independent people with 33 end-period samples. Fish oil is a mixture; no separate EPA/DHA causal estimates. Fenofibrate is a drug comparator. Not validation of an acute 6h response."),
 input_files=c(if(raw_mode)file.path(raw,"GSE27385_RAW.tar") else file.path(raw,"GSE27385_series_matrix.txt.gz"),file.path(raw,"GSE27385_family.soft.gz")),
 qc_limit=if(raw_mode)"Independent CEL RMA and PM intensity QC; ambiguous probe-to-gene mappings excluded. Independent RNA integrity and PBMC cell fractions not available." else "Published processed intensity matrix, not independent CEL QC. Ambiguous probe-to-gene mappings explicitly excluded.")
