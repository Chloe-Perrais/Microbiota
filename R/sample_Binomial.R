#' @title Sampling function for the abundances of species in different habitats
#' @param Nhabitats prob nindiviual_per_habitat
#' @description Functions that randomly sample abundance data in N habitat and return the abundance table
#'
#' @return dataframe
#' @export
sample_Binomial <- function(Nhabitats, prob, nindiviual_per_habitat){
  data <- rbinom(n=Nhabitats, p=prob, size=nindiviual_per_habitat)
  data <- rbind(data,nindiviual_per_habitat-data)
  colnames(data)=paste0("sample_",1:100)
  rownames(data)= c("OTU1","OTU2")
  return(data)
}
#' @examples sample_Binomial(Nhabitats=100, prob=0.5, nindiviual_per_habitat=5)