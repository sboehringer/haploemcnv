% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Geno2Diplo.R
\name{Geno2Haplo}
\alias{Geno2Haplo}
\title{Convert a genotype list into a diplotype list}
\usage{
Geno2Haplo(genotypes)
}
\arguments{
\item{genotypes}{A list. Needs to be in the format: GeneA1+GeneA2&GeneB1-GeneB2.}
}
\value{
A list with all diplotypes that are compatible with the individual's genotype. The output of a non-ambiguous case is as: GeneA1-GeneB1+GeneA2-GeneB2 and GeneA1-GeneB2+GeneA2-GeneB1.
}
\description{
To convert the genotype information into diplotype information.
}
