#' Loads GTEx metadata
#' 
#' @param SubjectPhenotypes path to SubjectPhenotypes file.
#' @param SubjectSampleMapping path to SubjectSampleMapping file.
#' @param SampleAttributes path to SampleAttributes file.
#' @return data.frame with merged metadata
#' @export

load_metadata <- function(
  SubjectPhenotypes,
  SubjectSampleMapping,
  SampleAttributes
  ) {
  SubjectPhenotypes <- read.delim(SubjectPhenotypes, fill=F)
  SubjectSampleMapping <- read.delim(SubjectSampleMapping, fill=F)
  SampleAttributes <- read.delim(SampleAttributes, fill=F)
  metadata <- merge(SubjectPhenotypes,SubjectSampleMapping, by='SUBJID')
  metadata <- merge(metadata,SampleAttributes, by='SAMPID')
  metadata
}

#' Generates a boxplot with the distribution of the values across the GTEx tissues
#' 
#' @param values data.frame with the values. rownames need to be the feature identifier (gene_id, exon_id,...) and colnames the sample ids.
#' @param metadata object returned by load_metadata()
#' @param gene.rpkm path to object gene.rpkm; load('/users/rg/freverter/GTExAnalysis/gtex_gc19_2015_01_12/rna_seq/gene_rpkm.RData')
#' @param gene_ids list of genes
#' @param log10_scale if TRUE (default), values are transfomed by log10
#' @param outlier.size size of outlier points in the boxplot
#' @return ggplot object
#' @export
#' @seealso \code{\link{load_metadata}}


ggtex_boxplot_tissues <- function(
  values,
  metadata,
  gene.rpkm,
  gene_ids,
  log10_scale = T,
  outlier.size = 1
  ) {
  if (values){
    
  }
  # check if gene.rpkm was set
  else if (gene.rpkm){
    if (!gene_ids){ # gene.rpkm requires a list of gene ids (ENSG) in order to get a subset of genes
      stop("gene_ids was not specified. If you provide the gene.rpkm, a list of genes must be specified using the parameter gene_ids",
           call. = FALSE)
    }
    # add a column with the ENSG without the suffix (e.g. ENSG00000223972.4 become ENSG00000223972)
    names_short <- as.character(lapply(strsplit(as.character(gene.rpkm$Name), split="\\."), "[",1))
    # do the same operation with the gene_ids list
    gene_ids_short <- as.character(lapply(strsplit(as.character(gene_ids), split="\\."), "[",1))
    # generate the subset of the genes
    genes <- names_short %in% gene_ids_short
    if (length(genes)==0) stop("None of the gene_ids was found in gene.rpkm$Name. Please, check if they have genes in common")
    # get their values excluding the first two columns (Name and Description)
    values <- gene.rpkm[genes,-c(1,2)]
    # and include the Description as rownames. It will be displayed in the plot
    rownames(values) <- gene.rpkm[genes,]$Description
  }

  # load values
  names(values) <- gsub("\\.","-", names(values))
  df <- melt(t(values))
  names(df) <- c('SAMPID','feature','value')
  df <- merge(df,metadata, by='SAMPID')
  
  # plot
  if (log10_scale) {
    p <- ggplot(df, aes_string('SMTSD', 'log10(value)', fill='SMTS'))
  } else { 
    p <- ggplot(df, aes_string('SMTSD', 'value', fill='SMTS'))
  }
  
  p <- p + geom_boxplot( outlier.size = outlier.size) +
    facet_grid(feature~SMTS, space="free", scales="free") +
    theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1, size=12),
          strip.text.x = element_text(angle = 90, size = 10, hjust = 0, vjust = 0.5),
          strip.text.y = element_text(angle = 0, size = 10, hjust = 0, vjust = 0.5)) +
    scale_fill_discrete(name="Tissue")
  
  p
}
  