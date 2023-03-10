% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Diplo2Geno.R
\name{Diplo2Geno}
\alias{Diplo2Geno}
\title{Convert a diplotype list into a genotype list}
\usage{
Diplo2Geno(diplotypes)
}
\arguments{
\item{diplotypes}{A list. Needs to be in the format: GeneA1-GeneB1+GeneA2-GeneB2, as obtained from \code{\link{combining_genes}}}
}
\value{
A list with all genotypes that are compatible with the individual's diplotype. The output of a non-ambiguous case is as: GeneA1+GeneA2&GeneB1+GeneB2
}
\description{
To convert the diplotype information into genotype information.
}
