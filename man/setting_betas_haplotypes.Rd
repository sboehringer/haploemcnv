% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setting_betas_haplotypes.R
\name{setting_betas_haplotypes}
\alias{setting_betas_haplotypes}
\title{Setting betas for all haplotypes based individual allele effect sizes}
\usage{
setting_betas_haplotypes(
  haplotypes,
  beta_alleles,
  gene_indication = TRUE,
  add_gene_indication = FALSE,
  combinations = NULL
)
}
\arguments{
\item{haplotypes}{A vector. Haplotypes and their frequencies.}

\item{beta_alleles}{An list. Overview of the effect sizes of each allele, each list element contains the alleles of a different gene. Can be obtained via \code{KIR_simulation_preparation}.}

\item{gene_indication}{A logical scalar. Whether or not it is indicated from which gene combination the allele or haplotype is (default is \emph{TRUE}).}

\item{add_gene_indication}{A logical scalar. Whether or not the indication from which gene combination the allele or haplotype is should be included (default is \emph{TRUE}).}

\item{combinations}{An optional list. Which allele combinations have full effect sizes, and what the loss in effect is of alleles not in this combination. Allows for the addition of haplotype effects. # Note that if \code{adj_scale} is 0 this list is ignored.}
}
\value{
A vector with the estimates beta for each haplotype combination
}
\description{
Setting betas for all haplotypes based individual allele effect sizes
}
