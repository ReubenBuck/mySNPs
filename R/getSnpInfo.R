#' Collect the genotype and allele frequency for a given snp
#' 
#' @param snp A vector of SNPs in AA/AB/BB format
#' @param pop A vector describing the population structure of the SNPs
#' @examples
#' snp1 <- rep(c("AA", "AB", "BB"), times = c(5,7,9))
#' pop <- sample(x = c("Treat", "Control"), size = 21, replace = TRUE)
#' getSnpInfo(snp1,pop)
#' @return A list with components \code{genotypes} and \code{allelFreqs}, which contain a breakdown of observed
#' genotypes and allele frequencies for each population
#' @export
getSnpInfo <- function(snp, pop = NULL){
  if (is.null(pop)){
    genotypes <- table(snp)
    freqs <- (genotypes[1] + 0.5*genotypes[2]) / sum(genotypes)
  }else{
    stopifnot(length(pop) == length(snp))
    genotypes <- table(pop, snp)
    freqs <- apply(genotypes, 1,
                   function(x){ (x[1] + .5*x[2]) / sum(x) })
  }
  list(genotypes = genotypes, alleleFreqs = freqs)
}










