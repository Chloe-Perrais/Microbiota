#' @title Compute variance
#' @param data_group abundance table with the number of reads per OTU per sample
#' @param last_position last letter of the sample
#' @description Function that give the variance of a group
#'
#' @return value
#' @export
Compute_var <- function(taxon, data_group, last_position){
  
  dataGRP_sum_final <- dplyr::filter(data_group, family==taxon) %>%
    pivot_longer(cols = !family, names_to = "Replicate", values_to = "Count")

  dataGRP_sum_final$Sample <- substr(dataGRP_sum_final$Replicate,1,last_position)
  dataGRP_sum_final$PCR <- as.factor(substr(dataGRP_sum_final$Replicate,last_position+2,last_position+2))
  dataGRP_sum_final$Extraction <- as.factor(ifelse(dataGRP_sum_final$PCR=="C",2,1))
 
#Linear model
  m0 <- lme4::lmer(log(Count+1)~1+(1|PCR)+(1|Extraction), data = dataGRP_sum_final)
  summary(m0)
## Extract variance components
  varcomp <- data.frame(lme4::VarCorr(m0))$vcov
  names(varcomp) <- data.frame(lme4::VarCorr(m0))$grp
  varcomp <- varcomp/sum(varcomp)
## Compute repeatability
#repeatability = varcomp[1]/sum(varcomp)
  return(varcomp)
}
#' @examples Compute_var(data_group=data.frame(x1_A=c(5:9),x1_B=c(5:9),x1_C=c(5:9),x2_A=c(15:19),x2_B=c(15:19),x2_C=c(15:19),family=c("ABC","ABC","ABC","ABC","ABC")), last_position=2)
