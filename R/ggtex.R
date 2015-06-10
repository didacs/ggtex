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