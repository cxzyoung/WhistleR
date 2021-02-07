#' @title Extract sequence-derived predictive features from interval-based data
#' @description A function to extract sequence features from the input \code{\link{GRanges}} object and the \code{\link{BSgenome}} object.
#' 
#' @param x A \link{GRanges} object for the genomic ranges to be annotated, the \code{width} of x must all be equal.
#' @param sequence A \link{BSgenome} or \link{XStringSet} object for the genome sequence. 
#' @param encoding Can be one of the following:
#' \describe{
#' \item{\strong{onehot}}{From the 5' end to the 3'end of x, 
#' each nucleotide position is coded by 4 indicators/dummy variables, 
#' where each dummy variable indicates that the position is equal to the base "A", "T", "C", and "G", respectively.}
#'
#' \item{\strong{iRNA}}{
#' Each nucleotide position is encoded by 4 variables, 
#' the first variable indicates that the nucleotide is purine (A or G), 
#' the second variable indicates that the nucleotide has an amino group (A or C), 
#' and the third variable indicates the formation of weak hydrogen bond (A or T), 
#' the fourth variable is calculated by the cumulative frequency of nucleotides from 
#' the leftmost position to the current position.
#' }
#' }
#' 
#' @details The function first extract DNA sequence within the genomic ranges defined by \code{x}. 
#' Then, the DNA strings are processed according to the selected sequence encoding method.
#' 
#' @return A \code{data.frame} object whose number of rows is the length of \code{x}, and the number of columns is 4 times the width of \code{x}.
#' The column types in the \code{data.frame} are all numeric.
#' 
#' @examples 
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' 
#' ## Define the Granges to be annotated:
#' set.seed(01)
#' 
#' X <- GRanges(rep(c("chr1", "chr2"), c(15, 15)),
#'              IRanges(c(sample(11874:12127, 15), sample(38814:41527,15)), width=5),
#'              strand=Rle(c("+", "-"), c(15, 15)))
#'    
#' ## Extract onehot encoded sequence features
#' seq_onehot <- sequenceDerivedFeatures(X, Hsapiens, encoding = "onehot")       
#' str(seq_onehot)                
#' 
#' ## Extract iRNA encoded sequence features
#' seq_iRNA <- sequenceDerivedFeatures(X, Hsapiens, encoding = "iRNA")
#' str(seq_iRNA)
#' 
#' @seealso 
#' 
#' \itemize{
#' \item{}{The \link{genomeDerivedFeatures} for extraction of genome-derived features.}
#' \item{}{The \link{topologyOnTranscripts} for calculation of the meta-tx topologies on transcripts of genes.}
#' }
#' 
#' @importFrom matrixStats rowCumsums
#' @importFrom IRanges Views
#' @importFrom Biostrings DNAStringSet
#' @export
sequenceDerivedFeatures <- function(x,
                             sequence,
                             encoding = c("onehot","iRNA")){
  encoding = match.arg(encoding)
  stopifnot(all(width(x) == width(x)[1]))
  
  ##Fill unknown attributes of x with the attributes of BSgenome
  if(anyNA(isCircular(x))){
    isCircularGenome <- isCircular(sequence)
    isCircular(x) <- isCircularGenome[names(isCircular(x))]
    rm(isCircularGenome)
  }
  if(anyNA(seqlengths(x))){
    seqLenGenome <- seqlengths(sequence)
    seqlengths(x) <- seqLenGenome[names(seqlengths(x))]
    rm(seqLenGenome)
  }
  
  N = width(x)[1]
  sequences <- as.character(DNAStringSet(Views(sequence, x)))
  sequence_M <- matrix( unlist( strsplit(sequences,"") ), ncol =  N, byrow = TRUE)
  rm(sequences)
  
  if(encoding == "onehot"){
  seqfeatures <- onehot_encode(sequence_M)
  }else if(encoding == "iRNA"){
  seqfeatures <- iRNA_encode(sequence_M)
  }
  seqfeatures <- as.data.frame(seqfeatures)
  return(seqfeatures)
}

onehot_encode <- function(sequence_M){
  A_M <- sequence_M == "A"
  colnames(A_M) <- paste0("A_", seq_len(ncol(A_M)))
  T_M <- sequence_M == "T"
  colnames(T_M) <- paste0("T_", seq_len(ncol(T_M)))
  C_M <- sequence_M == "C"
  colnames(C_M) <- paste0("C_", seq_len(ncol(C_M)))
  G_M <- sequence_M == "G"
  colnames(G_M) <- paste0("G_", seq_len(ncol(G_M)))
  N <- ncol(sequence_M)
  rm(sequence_M)
  onehotfeatures <- cbind(A_M, T_M, C_M, G_M) 
  rm(A_M, T_M, C_M, G_M)
  onehotfeatures <- as.data.frame(onehotfeatures[,rep(seq_len(N), each=4) + rep(c(0,N,2*N,3*N),N)])
  onehotfeatures <- apply(onehotfeatures, 2, as.numeric)
  return(onehotfeatures)
}

iRNA_encode <- function(sequence_M){
  N <- ncol(sequence_M)
  purine_M <- sequence_M == "A" | sequence_M == "G" 
  colnames(purine_M) <-  paste0("purine_",seq_len(N))
  amino_M <- sequence_M == "A" | sequence_M == "C" 
  colnames(amino_M) <-  paste0("aminoGroup_",seq_len(N))
  weakHyb_M <- sequence_M == "A" | sequence_M == "T" 
  colnames(weakHyb_M) <- paste0("weakHydrogenBonds_",seq_len(N))
  
  cumFreq_A <- rowCumsums(matrix(as.numeric( sequence_M == "A" ), ncol = N, byrow = FALSE))
  cumFreq_T <- rowCumsums(matrix(as.numeric( sequence_M == "T" ), ncol = N, byrow = FALSE))
  cumFreq_C <- rowCumsums(matrix(as.numeric( sequence_M == "C" ), ncol = N, byrow = FALSE))
  cumFreq_G <- rowCumsums(matrix(as.numeric( sequence_M == "G" ), ncol = N, byrow = FALSE))
  
  cumFreq_combined <- matrix(0, ncol = N, nrow = nrow(sequence_M))
  cumFreq_combined[sequence_M == "A"] <- cumFreq_A[sequence_M == "A"]
  cumFreq_combined[sequence_M == "T"] <- cumFreq_T[sequence_M == "T"]
  cumFreq_combined[sequence_M == "C"] <- cumFreq_C[sequence_M == "C"]
  cumFreq_combined[sequence_M == "G"] <- cumFreq_G[sequence_M == "G"]
  
  cumFreq_combined <- t(t(cumFreq_combined) / seq_len(N))
  colnames(cumFreq_combined) <-  paste0("cumulativeFrequency_",seq_len(N))
  seqfeatures <- as.data.frame(cbind(cbind(purine_M,amino_M,weakHyb_M),as.data.frame(cumFreq_combined)))
  seqfeatures <- apply(seqfeatures, 2, as.numeric)
  return(seqfeatures)
}

## To Do: add more sequence feature encoding methods.

