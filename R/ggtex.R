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
#' @param gene.rpkm gene.rpkm object; load('/users/rg/freverter/GTExAnalysis/gtex_gc19_2015_01_12/rna_seq/gene_rpkm.RData')
#' @param transcript.rpkm transcript.rpkm object; load('/users/rg/freverter/GTExAnalysis/gtex_gc19_2015_01_12/rna_seq/transcript_rpkm.RData')
#' @param gene_ids list of genes
#' @param transcript_ids list of genes
#' @param outlier.size size of outlier points in the boxplot
#' @return ggplot object
#' @export
#' @seealso \code{\link{load_metadata}}


ggtex_boxplot_tissues <- function(
  values=F,
  metadata,
  gene.rpkm=F,
  transcript.rpkm=F,
  gene_ids=F,
  transcript_ids=F,
#  log10_scale = T,
  outlier.size = 1
  ) {
  
  # values by default, they will be redefined if transcript.rpkm is in input
  fill='SMTS'
  name="Tissue"
  
  if (class(values)!="data.frame") {
    # values is going to be obtained from the loaded object
    
    # check input parameters
    if ((class(gene.rpkm)!="data.frame") & (class(transcript.rpkm)!="data.frame")){
      stop("no values were provided, please, use values, gene.rpkm or transcript_ids.")
    }
    if ((class(gene.rpkm)=="data.frame") & (class(transcript.rpkm)=="data.frame")){
      stop("gene.rpkm or transcript.rpkm, but not both please.")
    }
    if ((class(gene_ids)=="logical") & (class(transcript_ids)=="logical")) {
      stop("no genes or transcripts list, please, use gene_ids or transcript_ids")
    }
    
    # gene.rpkm as input
    if (class(gene.rpkm)=="data.frame"){
      if (class(gene_ids)=="logical"){ # gene.rpkm requires a list of gene ids (ENSG) in order to get a subset of genes
        stop("gene_ids was not specified. Since you provided the gene.rpkm, a list of genes must be specified with the parameter gene_ids.")
      }
      # ENSG without the suffix (e.g. ENSG00000223972.4 become ENSG00000223972)
      names_short <- as.character(lapply(strsplit(as.character(gene.rpkm$Name), split="\\."), "[",1))
      # same operation with the gene_ids list
      gene_ids_short <- as.character(lapply(strsplit(as.character(gene_ids), split="\\."), "[",1))
      # generate the subset of the genes
      genes <- names_short %in% gene_ids_short
      # check if genes is not empty  
      if (length(genes)==0) stop("None of the gene_ids was found in gene.rpkm$Name. Please, check if they have any genes in common")
      # get their values excluding the first two columns (Name and Description)
      values <- gene.rpkm[genes,-c(1,2)]
      # and include the Description as rownames. It will be displayed in the plot
      rownames(values) <- gene.rpkm[genes,]$Description
    } 
    
    # transcript.rpkm in input
    else if (class(transcript.rpkm)=="data.frame"){
      column = F
      ids_list = F
      if (class(gene_ids)!="logical"){
        column <- 'Gene_Symbol'
        ids_list <- gene_ids
      } else if (class(transcript_ids)!="logical"){ 
        column = 'TargetID'
        ids_list <- transcript_ids
      }
      if (class(column)=="logical") stop("no genes or transcripts list, please, use gene_ids or transcript_ids")
      if (class(ids_list)=="logical") stop("no genes or transcripts list, please, use gene_ids or transcript_ids")
      
      # ids without the suffix (e.g. ENSG00000223972.4 become ENSG00000223972)
      object_ids_short <- as.character(lapply(strsplit(as.character(transcript.rpkm[,column]), split="\\."), "[",1))
      # same operation with the gene_ids list
      input_ids_short <- as.character(lapply(strsplit(as.character(ids_list), split="\\."), "[",1))
      # generate the subset of the genes
      features <- object_ids_short %in% input_ids_short
      # check if genes is not empty  
      if (length(features)==0) stop("None of the ids in input was found in the values object. Please, check if they have any id in common")
      values <- transcript.rpkm[features,-c(1:4)]
      # and include the Description as rownames. It will be displayed in the plot
      rownames(values) <- transcript.rpkm[features,'TargetID']
      
    }
  }
  # transform values
  # replace "." by "-" in colnames
  names(values) <- gsub("\\.","-", names(values))
  # data.frame of values in long format
  df <- melt(t(values))
  names(df) <- c('SAMPID','feature','value')
  df <- merge(df, metadata, by='SAMPID')
  
  if (class(transcript.rpkm)=="data.frame") {
    rows <- match( df$feature, transcript.rpkm$TargetID )
    df$Gene_Symbol <- factor( transcript.rpkm[rows, 'Gene_Symbol'] )
    df$feature <- reorder(df$feature, order(df$Gene_Symbol))
    fill <- 'Gene_Symbol'
    name <- 'Gene'
  }
  
  # plot
  p <- ggplot(df, aes_string('SMTSD', 'value', fill=fill))
  p <- p + geom_boxplot( outlier.size = outlier.size) +
    facet_grid(feature~SMTS, space="free_x", scales="free") +
    theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1, size=12),
          strip.text.x = element_text(angle = 90, size = 10, hjust = 0, vjust = 0.5),
          strip.text.y = element_text(angle = 0, size = 10, hjust = 0, vjust = 0.5)) +
    scale_fill_discrete(name=name)
  
  p
}


#' Generates a boxplot with the distribution of the values across the GTEx tissues, including description of features
#' 
#' @param values data.frame with the values. rownames need to be the feature identifier (gene_id, exon_id,...) and colnames the sample ids.
#' @param metadata object returned by load_metadata()
#' @param exons exons
#' @param junctions junctions
#' @export
#' @seealso \code{\link{load_metadata}}


boxplot_features <- function(
  values,
  metadata,
  exons,
  junctions){
  
  # transform values
  # replace "." by "-" in colnames
  names(values) <- gsub("\\.","-", names(values))
  # data.frame of values in long format
  df <- melt(t(values))
  names(df) <- c('SAMPID','feature','value')
  df <- merge(df, metadata, by='SAMPID')
  
  p1 <- ggplot(df, aes_string('SMTSD', 'value', fill='SMTS'))
  p1 <- p1 + geom_boxplot( outlier.size = 1) +
    facet_grid(feature~SMTS, space="free_x", scales="free") +
    theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1, size=12),
          strip.text.x = element_text(angle = 90, size = 10, hjust = 0, vjust = 0.5),
          strip.text.y = element_text(angle = 0, size = 10, hjust = 0, vjust = 0.5)) +
    scale_fill_discrete(name='Tissues')
  
#  j <- read.table(junctions, header=T)
  j <- within(j, {
    chr <- as.numeric(lapply(strsplit(as.character(j$TargetID), split = '_'),"[",1))
    start <- as.numeric(lapply(strsplit(as.character(j$TargetID), split = '_'),"[",2))
    end <- as.numeric(lapply(strsplit(as.character(j$TargetID), split = '_'),"[",3))
    midpoint <- start + (end - start) / 2
  })
 # e <- read.table(exons, header=F, col.names = c('ID','chr','start','end','strand'))
  
  p2 <- ggplot()
  p2 <- p2 + geom_linerange(data=e[-1,], aes(x=0, ymin=start, ymax=end), size=10, color = "grey50")
  p2 <- p2 + geom_segment(data=transform(j, TargetID=NULL), aes(x=0, xend=1, y=start, yend=midpoint), size=.5)
  p2 <- p2 + geom_segment(data=transform(j, TargetID=NULL), aes(x=1, xend=0, y=midpoint, yend=end), size=.5)
  p2 <- p2 + geom_segment(data=j, aes(x=0, xend=1, y=start, yend=midpoint), size=1.5, color = "red")
  p2 <- p2 + geom_segment(data=j, aes(x=1, xend=0, y=midpoint, yend=end), size=1.5, color = "red")
  p2 <- p2 + coord_flip() + xlim(-2,3) + facet_grid(TargetID~.)
  
  grid.arrange(p1, p2, ncol=2, main = "Main title", widths=unit.c(unit(0.65, "npc"), unit(0.35, "npc")))
  
  
}