% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KIR2DS4_combining.R
\name{KIR2DS4_combining}
\alias{KIR2DS4_combining}
\title{Combining the KIR2DS4 and KIR2DS4N KIR genes into one gene}
\usage{
KIR2DS4_combining(KIR2DS4, KIR2DS4N)
}
\arguments{
\item{KIR2DS4}{A vector. The genotyping data of KIR2DS4.}

\item{KIR2DS4N}{A vector. The genotyping data of KIR2DS4N.}
}
\value{
A vector with the combined genotyping result of the two genes. The combination made for each individual is also printed in the Rstudio console.
}
\description{
To combine the two KIR2DS4 genes into one. The gene copies of the two genes are simply added together.
}
\details{
If any of the two genes has a `POS` genotype call, than the combined gene also gets the `POS` genotype.
}
