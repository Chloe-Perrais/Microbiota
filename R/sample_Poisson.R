#' @title Sampling function for Poisson data
#' @param sample count lambda
#' @description Functions that randomly samples Poisson data and returns a dataframe
#'
#' @return dataframe
#' @export
sample_Poisson <- function(sample, count, lambda){
  data <- data.frame(sample=rep(sample), count=rpois(n=count, lambda = lambda))
  return(data)
}
#' @examples sample_Poisson(sample=c("sample_1","sample_2"), count=1000, lambda= c(10,100))