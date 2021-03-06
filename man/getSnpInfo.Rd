% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSnpInfo.R
\name{getSnpInfo}
\alias{getSnpInfo}
\title{Collect the genotype and allele frequency for a given snp}
\usage{
getSnpInfo(snp, pop = NULL)
}
\arguments{
\item{snp}{A vector of SNPs in AA/AB/BB format}

\item{pop}{A vector describing the population structure of the SNPs}
}
\value{
A list with components \code{genotypes} and \code{allelFreqs}, which contain a breakdown of observed
genotypes and allele frequencies for each population
}
\description{
Collect the genotype and allele frequency for a given snp
}
\examples{
snp1 <- rep(c("AA", "AB", "BB"), times = c(5,7,9))
pop <- sample(x = c("Treat", "Control"), size = 21, replace = TRUE)
getSnpInfo(snp1,pop)
}

