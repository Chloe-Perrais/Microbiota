#' @title Find chimeras
#' @param fordada abundance table
#' @description Functions that finds and cleans chimeras
#'
#' @export
find_chimeras <- function(fordada, name){
  print("find_chimeras")
  print(colnames(fordada))
  # Remove chimeras
  ## Tune data structure
  samples <- colnames(fordada)[seq(1+which(colnames(fordada)=="observation_sum"), ncol(fordada))]
  
  print("samples:")
  print(samples)
  fordada.m <- data.matrix(fordada[ , seq(1+which(colnames(fordada)=="observation_sum"), ncol(fordada))])
  rownames(fordada.m) <- fordada$seed_sequence
  fordada.t <- t(fordada.m)

  ## Remove chimeras with isBimeraDenovo
  bimeras.v <- dada2::isBimeraDenovo(dada2::getUniques(fordada.t), minFoldParentOverAbundance=4, verbose=TRUE)
  bimeras.df <- data.frame(x=as.logical(bimeras.v), seq=names(bimeras.v))
  head(bimeras.df)

  ## Export sequences that are not chimeras
  write.table(fordada[bimeras.df$x==FALSE,], file = file.path(here::here(), "data", "derived_data", paste(name,"_cleaned_abundance.tsv",sep="")), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)#writing of filtered abundance table (without identified chimeras)
  
  ## Export sequences that are chimeras
  write.table(fordada[bimeras.df$x==TRUE,], file = file.path(here::here(), "data", "derived_data", paste(name,"_chimeras_list.tsv", sep="")), sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)#downloading of identified chimeras
  
  ## Summary file
  number_of_total_clusters = nrow(fordada)
  number_of_bimeras = nrow(fordada[bimeras.df$x==TRUE,])
  number_of_nobimeras = nrow(fordada[bimeras.df$x==FALSE,])
  summary.df = data.frame(number_of_total_clusters, number_of_bimeras, number_of_nobimeras)
  write.table(summary.df,file=file.path(here::here(), "data", "derived_data", paste(name,"_chimeras_summary.txt",sep="")),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  
}
#' @examples find_chimeras(fordada=fordada, name=name)