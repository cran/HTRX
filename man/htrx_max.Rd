% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/htrx_max.R
\name{htrx_max}
\alias{htrx_max}
\alias{htr_max}
\title{Maximum independent features}
\usage{
htrx_max(nsnp, n_haps = 3^nsnp - 1, cap = 40)

htr_max(nsnp, cap = 40)
}
\arguments{
\item{nsnp}{a positive integer giving the number of single nucleotide polymorphisms (SNPs) included in the haplotypes.}

\item{n_haps}{a positive integer giving the number of haplotypes,
which is also the number of columns of the HTRX or HTR matrix.
By default, n_haps=3^nsnp-1}

\item{cap}{a positive integer which manually sets the maximum number of independent features.
By default, cap=40.}
}
\value{
a positive integer giving the maximum number of independent features to be included in the analysis.
}
\description{
The maximum number of independent features in principle from
all the haplotypes generated by Haplotype Trend Regression with eXtra flexibility (HTRX) (\code{\link{htrx_max}}) or Haplotype Trend Regression (HTR) (\code{\link{htr_max}}).
}
\details{
The maximum number of independent features in principle is \ifelse{html}{\out{2<sup>nsnp</sup>}}{\eqn{2^nsnp}}-1.
You can also manually set the upper limit of independent features (by setting "cap") that can be included in the final HTRX or HTR model.
}
\examples{
## the maximum number of independent haplotypes consisted of 4 SNPs from HTRX
htrx_max(nsnp=4,n_haps=(3^4-1))

## the maximum number of independent haplotypes consisted of 6 SNPs from HTR
htr_max(nsnp=6)
}