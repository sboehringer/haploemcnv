% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD.R
\name{LD}
\alias{LD}
\title{linkage disequilibrium calculation}
\usage{
LD(EM_out = NULL, haplotypes = NULL, three_plus = FALSE)
}
\arguments{
\item{EM_out}{An EM-algorithm analysis of at least two genes, as provided by \code{\link{EM_algorithm}}.}

\item{haplotypes}{An optional vector. Haplotypes and their probabilities, must be haplotypes consisting of at least two loci.}

\item{three_plus}{A logical scalar. Whether or not the linkage disequilibrium is calculated between the K-1 genes and gene K, or between gene 1 and gene 2 (\code{FALSE} is default).}
}
\value{
An matrix with the D', delta_{ij}, the R2, the observed and the expected haplotype frequency as separate columns and all possible gene combinations as different rows.
}
\description{
Linkage disequilibrium is calculated between two genes or between a combination of genes and the latest added gene.
}
