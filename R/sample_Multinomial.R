#' @title Sampling function for the abundances of species in different habitats
#' @param Nhabitats number of habitats (could be either a vector of values or a single value that will be recycle)
#' @param prob frequency of all OTUs (could be either a vector of values or a single value that will be recycle)
#' @param mean_nindiviuals_per_habitat average number of individuals per habitat (could be either a vector of values or a single value that will be recycle)
#' @description Functions that randomly sample abundance data in N habitat and return the abundance table
#'
#' @return dataframe
#' @export
sample_Multinomial <- function(Nhabitats, prob, mean_nindiviuals_per_habitat){
  if(length(mean_nindiviuals_per_habitat)==1){
    count <- rpois(n=Nhabitats, lambda = mean_nindiviuals_per_habitat)
    data <- t(mc2d::rmultinomial(n=Nhabitats, p=prob, size=count))

  }else{
    data <- t(mc2d::rmultinomial(n=Nhabitats, p=prob, size=mean_nindiviuals_per_habitat))
  }
  colnames(data)=paste0("sample_", 1:Nhabitats)
  data <- data.frame(otu=paste0("OTU", 1:length(prob)), data)
  #rownames(data)= paste0("OTU", 1:length(prob))
  return(data)
}
#' @examples sample_Multinomial(Nhabitats=100, prob=0.5, mean_nindiviuals_per_habitat=5)