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
#' @param metadata object obtained by load_metadata()
#' @param log10_scale if TRUE (default), values are transfomed by log10
#' @param outlier.size size of outlier points in the boxplot
#' @return ggplot object
#' @export
#' @seealso \code{\link{load_metadata}}


ggtex_boxplot_tissues <- function(
  values,
  metadata,
  log10_scale = T,
  outlier.size = 1
  ) {
  # load values
  names(values) <- gsub("\\.","-", names(values))
  df <- melt(t(values))
  names(df) <- c('SAMPID','feature','value')
  df <- merge(df,metadata, by='SAMPID')
  
  # plot
  
  if (log10_scale) {
    p <- ggplot(df, aes(SMTSD, log10(value), fill=SMTS))
  } else { 
    p <- ggplot(df, aes(SMTSD, value, fill=SMTS))
  }
  
  p + geom_boxplot( outlier.size = outlier.size) +
    facet_grid(feature~SMTS, space="free", scales="free") +
    theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1, size=12),
          strip.text.x = element_text(angle = 90, size = 10, hjust = 0, vjust = 0.5),
          strip.text.y = element_text(angle = 90, size = 10, hjust = 0, vjust = 0.5)) +
    scale_fill_discrete(name="Tissue")
  
  p
}
  