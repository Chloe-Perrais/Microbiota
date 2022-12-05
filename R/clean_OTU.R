#' @title Clean the abundance data
#' @param abundance_table abundance table
#' @param affiliation_table affiliation table
#' @param thresh % threshold below which 'unknown' are ignored
#' @param name name will be used to import files and write names before extension of the different output tables
#' @param removecol remove columns containing this specific word
#' @description Work on raw abundance and multihits text files from Frogs.
#'
#'Correct the taxonomy of each cluster of the abundance tsv file tagged 'multi-affiliation', from decision rules applied on the cleaned multihits tsv file:
#'  
#'  Clean the taxonomy (7 ranks: Kingdom,Phylum,Class,Order,Family,Genus,Species) of each hit of the two files for:
#'   - the special characters '[]()_"0123456789.' in all 7 taxonomic ranks.
#' - 'mitochondria' (in the Family rank) and 'chloroplast' (in the Class rank), which are replaced by plantae in all 7 taxonomic ranks.
#' - 'unidentified', 'sp.' and 'bacterium', which are replaced by 'unknown species' in the Species rank .
#' - unexpected number of words in all 7 taxonomic ranks.
#' - for the species and genus ranks only: if the number of 'unknown species' hits < % threshold then ignore them.
#' - for all taxonomic ranks: 
#' - if there is a unique taxonomic name, then rename by this taxonomic name.
#' - if there is a unique taxonomic name and an 'unknown', then rename by 'unknown'.
#' - if there are two or more taxonomic names, then name 'multi-affiliation'.
#' Identify chimeric sequences using the isBimeraDenovo(...) function in the dada2 library
#' Export an abundance table tsv file cleaned for taxonomy and filtered for chimeras.
#'
#' @export
clean_OTU <- function(abundance_table = "Galaxy21-[FROGS_BIOM_to_TSV__abundance.tsv].tsv", affiliation_table = "Galaxy22-[FROGS_BIOM_to_TSV__multi-affiliations.tsv].tsv",   thresh = 0.5, name = "2022JRLAT-16S", removecol=NULL){
  
  ## Import multi-affiliation data
  mt <- data.table::fread(file=file.path(here::here(), "data", "raw_data", affiliation_table), sep="\t", header=TRUE) # to modify (name of the multihits tsv file)
  dim(mt)
  names(mt)
  mh <- data.table::data.table(bidouille=mt$blast_taxonomy, blast_taxonomy=mt$blast_taxonomy, observation_name=mt$`#observation_name`)
  names(mh)
  ## Remove file
  rm(mt)
  
  ## Import abundance table
  at <- data.table::fread(file=file.path(here::here(), "data", "raw_data", abundance_table),sep="\t",header=TRUE)# to modify (name of the abundance file)
  if(!is.null(removecol)){
    at <- at %>% 
      dplyr::select(-contains(removecol))
    print("Remaining colnames")
    print(colnames(at))
  }
    ab <- data.table::data.table(bidouille=at$`blast_taxonomy`, blast_taxonomy=at$`blast_taxonomy`, observation_name=at$observation_name)

  # Cleaning of multihits file
  ## Cleaning for the special characters []()_"0123456789.
  mh$bidouille <- gsub("[\\[\\]()_\"0-9\\.]", "", mh$bidouille, perl=TRUE)
  mh$bidouille <- gsub(" +;", ";", mh$bidouille, perl=TRUE)
  
  ## Deconstruction of the taxonomic column ('blast_taxonomy') in 7 taxonomic fields
  ## Warnings message when no data on taxonomy
  mh[mh$bidouille=="no data",]
  #print(names(mh))
  mh <- mh %>% tidyr::separate(col=bidouille, into=c("tx1","tx2", "tx3", "tx4", "tx5", "tx6", "tx7"), sep=";")
  #print(names(mh))

  #mh <- cbind(toto, mh[, c("observation_name" ,"blast_taxonomy")]) # observation name is the cluster ID and useful to link with the abundance tsv file
  #rm(toto)
  
  ## Replacement of 'chloroplast' and 'mitochondria' by 'plantae' in all 7 taxonomic fields
  ## Remove one "grepl("Chloroplast", mh$tx4, ignore.case=TRUE)"
  mh[grepl("Chloroplast", mh$tx4, ignore.case=TRUE)|grepl("Chloroplast", mh$tx4, ignore.case=TRUE)|grepl("Mitochondria", mh$tx5, ignore.case=TRUE), 1:7] <- c("Plantae")
  
  ## Overwriting of the species field with 'unknown species' if presence of 'incertae sedis' in any of upper taxonomic fields 
  ## Overwriting of the genus field with 'unknown genus'  if presence of 'incertae sedis' in any of upper taxonomic fields
  ## NB: no incertae sedis in JRL-AT dataset
  mh[!(grepl("incertae sedis", mh$tx2, ignore.case=TRUE)|grepl("incertae sedis", mh$tx3,ignore.case=TRUE)|grepl("incertae sedis", mh$tx4, ignore.case=TRUE)|grepl("incertae sedis", mh$tx5, ignore.case=TRUE)|grepl("incertae sedis", mh$tx6, ignore.case=TRUE))&grepl("incertae sedis", mh$tx7, ignore.case=TRUE), 7]<-c("unknown species")
  mh[grepl("incertae sedis", mh$tx2, ignore.case=TRUE)|grepl("incertae sedis", mh$tx3, ignore.case=TRUE)|grepl("incertae sedis", mh$tx4, ignore.case=TRUE)|grepl("incertae sedis", mh$tx5, ignore.case=TRUE)|grepl("incertae sedis", mh$tx6, ignore.case=TRUE), 6:7]<-list("unknown genus","unknown species")
  
  ## Cleaning for more than 1 word in the 6 first taxonomic fields (except if 'unknown' is found)
  ### Prior: removal of 'Candidatus' in each field
  mh <- as.data.frame(mh)
  for(j in 1:6) {
    mh[, j] <- gsub("^Candidatus ", "", mh[,j ])
    my.stuff <- as.data.frame(mh[, j])
    tf <- !grepl("unknown", my.stuff[, 1], ignore.case=TRUE)
    mh[which(tf),j] <- gsub(" +.*$", "", mh[which(tf), j], perl=TRUE)
  }
  mh <-  data.table::as.data.table(mh)
  
  ## Cleaning for more than 2 words in the 7th taxonomic field
  ## Prior: removal of 'Candidatus' in the species field
  mh$tx7 <- gsub("^Candidatus ","", mh$tx7, perl=TRUE)
  mh$tx7 <- gsub(" +$", "", gsub("^(([^ ]+ ?){1,2}).*$","\\1",mh$tx7,perl=TRUE),perl=TRUE)
  
  ## Replacement in the species field of 'unidentified' or 'metagenome' anywhere in the field by 'unknown species' 
  ## Replacement in the species field of 'sp' at the end of the field by 'unknown species'
  ## Replacement in the species field of 'bacterium' or 'uncultured' at the beginning of the field by 'unknown species'
  ## Replacement of the species field by 'unknown species' if empty
  mh$tx7 <- gsub("^.*(unidentified|metagenome).*$","unknown species",mh$tx7,perl=TRUE)
  mh$tx7 <- gsub("^.*(sp|eae bacterium|ales bacterium)$","unknown species",mh$tx7,perl=TRUE)
  mh$tx7 <- gsub("^(uncultured|bacterium).*$","unknown species", mh$tx7, perl=TRUE)
  mh$tx7[grepl("^ *$", mh$tx7)] <-"unknown species"
  
  ## Overwriting of the species field with 'unknown species' if it begins by a lower case 
  todo <- grepl("^[a-z]", mh$tx7) & !grepl("unknown species", mh$tx7)
  mh$tx7[which(todo)] <- "unknown species"
  
  ## Overwriting of the genus field  (if unknown) with the first word of the species field (if not unknown) if empty 
  todo <- mh$tx6=="unknown genus" & grepl(" ", mh$tx7) & !grepl("unknown species", mh$tx7)
  mh$tx6[which(todo)] <- gsub(" .*$", "", mh$tx7[which(todo)])
  
  ## Ignoring of 'unknown species' or 'unknown genus' under a threshold
  mh$clus <- as.integer(gsub("Cluster_", "", mh$observation_name))
  mh$blast_taxonomy <- NULL
  
  ## Add a column with the number of rows per cluster in affiliation table
  my.genus.result <- data.table::as.data.table(unique(mh %>% dplyr::add_count(clus, tx1, tx2, tx3, tx4, tx5, tx6, observation_name)))
  head(my.genus.result)
  my.genus.result$ntot <-  data.table::as.data.table(unique(mh %>% dplyr::add_count(observation_name)))$n
  head(my.genus.result)

  data.table::setkey(my.genus.result, clus, n)
  my.genus.result$prop.genus <- my.genus.result$n/my.genus.result$ntot
  my.genus.result$ignore.genus <-FALSE
  my.genus.result$ignore.genus <- my.genus.result$ignore.genus | (grepl("unknown", my.genus.result$tx6, perl=TRUE, ignore.case=TRUE) & my.genus.result$prop.genus < thresh)
  #nrow(my.genus.result)
  #names(my.genus.result)
  
  my.species.result <- data.table::as.data.table(unique(mh %>% dplyr::add_count(clus, tx1, tx2, tx3, tx4, tx5, tx6, tx7, observation_name)))
  my.species.result$ntot <- data.table::as.data.table(unique(mh %>% dplyr::add_count(observation_name)))$n
  data.table::setkey(my.species.result, clus, n)
  my.species.result$prop.species <- my.species.result$n/my.species.result$ntot
  my.species.result$ignore.species <- FALSE
  my.species.result$ignore.species <- my.species.result$ignore.species | (grepl("unknown", my.species.result$tx7, perl=TRUE,ignore.case=TRUE) & my.species.result$prop.species < thresh)
  #nrow(my.species.result)
  #names(my.species.result)
  
  ## Merge datasets together
  my.first.result <- merge(my.genus.result, my.species.result, by = names(my.species.result)[1:8])[, c(1:8,13,18)]
  
  my.second.result <- data.table::data.table(observation_name=unique(my.first.result$observation_name))
  for(j in 1:7) {
    tx <- as.character(as.data.frame(my.first.result)[,j])
    tmp <- apply(my.second.result, 1, function(r){
      v <- unique(tx[(j==6 | !my.first.result$ignore.genus) & my.first.result$observation_name==r[1]])
      v <- unique(tx[(j==7 | !my.first.result$ignore.species) & my.first.result$observation_name==r[1]])
      if (length(v)==1) {v[1]} 
      else if ((length(v)==2)) {
        unk <- grepl("unknown",v)
        unkV <- unique(v[unk])
        if (length(unkV) == 0) {
          "Multi-affiliation"
        } else {
          unkV[1]
        }
      } else {
        "Multi-affiliation"
      }
    })
    my.second.result <- cbind(my.second.result, tmp)
  }
  #nrow(my.second.result)
  
  ## Reconstruction of the new full taxonomy from each taxonomic field for each cluster
  my.mh.result <- data.table::data.table(observation_name=my.second.result$observation_name,  blast_taxonomy= apply(my.second.result, 1, function(r) {
    paste(r[2:8], sep=";", collapse=";")
  }))
  #nrow(my.mh.result)
  dim(my.mh.result)
  
  # Cleaning of abundance file
  ## Cleaning for some special characters
  ab$bidouille <- gsub("[\\[\\]()_\"0-9\\.]","", ab$bidouille, perl=TRUE)
  ab$bidouille <- gsub(" +;", ";", ab$bidouille, perl=TRUE)
 
  ## Deconstruction of the taxonomic column in 7 taxonomic fields
  ab <- ab %>% tidyr::separate(col=bidouille, into=c("tx1"," tx2", "tx3", "tx4", "tx5", "tx6", "tx7"), sep=";")
  #ab <- cbind(toto, ab[, c("observation_name" ,"blast_taxonomy")])
  #rm(toto)
  
  ## Replacement of 'chloroplast' and 'mitochondria' by 'plantae' in all 7 taxonomic fields
  ab[grepl("Chloroplast", ab$tx4,ignore.case=TRUE)|grepl("Mitochondria",ab$tx5,ignore.case=TRUE),1:7]<-c("Plantae")
  
  ## Overwriting of the species field with 'unknown species' if presence of 'incertae sedis' in any of upper taxonomic fields 
  ## Overwriting of the genus field with 'unknown genus'  if presence of 'incertae sedis' in any of upper taxonomic fields
  ab[!(grepl("incertae sedis", ab$tx2, ignore.case=TRUE)|grepl("incertae sedis", ab$tx3, ignore.case=TRUE)|grepl("incertae sedis", ab$tx4, ignore.case=TRUE)|grepl("incertae sedis", ab$tx5, ignore.case=TRUE)|grepl("incertae sedis", ab$tx6,ignore.case=TRUE))&grepl("incertae sedis", ab$tx7,ignore.case=TRUE), 7]<-c("unknown species")
  ab[grepl("incertae sedis", ab$tx2, ignore.case=TRUE)|grepl("incertae sedis",ab$tx3, ignore.case=TRUE)|grepl("incertae sedis", ab$tx4, ignore.case=TRUE)|grepl("incertae sedis", ab$tx5, ignore.case=TRUE)|grepl("incertae sedis", ab$tx6, ignore.case=TRUE), 6:7] <- list("unknown genus", "unknown species")
  
  ## Cleaning for more than 1 word in the 6 first taxonomic fields (except if 'unknown' is found)
  ### Prior: removal of 'Candidatus' in each field
  ab <- as.data.frame(ab)
  for(j in 1:6) {
    ab[,j] <- gsub("^Candidatus ", "", ab[,j])
    my.stuff <- as.data.frame(ab[,j])
    tf <- !grepl("unknown",my.stuff[,1],ignore.case=TRUE)
    ab[which(tf),j] <- gsub(" +.*$","", ab[which(tf),j], perl=TRUE)
  }
  ab <-  data.table::as.data.table(ab)
  
  ## Cleaning for more than 2 words in the 7th taxonomic field
  ## Prior: removal of 'Candidatus' in the species field
  ab$tx7<-gsub("^Candidatus ","",ab$tx7,perl=TRUE)
  ab$tx7<-gsub(" +$", "", gsub("^(([^ ]+ ?){1,2}).*$","\\1",ab$tx7,perl=TRUE),perl=TRUE)
  
  ## Replacement in the species field of 'unidentified' or 'metagenome' anywhere in the field by 'unknown species' 
  ## Replacement in the species field of 'sp' at the end of the field by 'unknown species'
  ## Replacement in the species field of 'bacterium' or 'uncultured' at the beginning of the field by 'unknown species'
  ## Replacement of the species field by 'unknown species' if empty
  ab$tx7<-gsub("^.*(unidentified|metagenome).*$","unknown species",ab$tx7,perl=TRUE)
  ab$tx7<-gsub("^.*(sp|eae bacterium|ales bacterium)$","unknown species",ab$tx7,perl=TRUE)
  ab$tx7<-gsub("^(uncultured|bacterium).*$","unknown species",ab$tx7,perl=TRUE)
  ab$tx7[grepl("^ *$",ab$tx7)]<-"unknown species"
  
  ## Overwriting of the species field with 'unknown species' if it begins by a lower case 
  todo <- grepl("^[a-z]", ab$tx7) & !grepl("unknown species", ab$tx7)
  ab$tx7[which(todo)] <- "unknown species"
  
  ## Overwriting of the genus field (if unknown) with the first word of the species field (if not unknown) if empty 
  todo <- ab$tx6=="unknown genus" & grepl(" ", ab$tx7) & !grepl("unknown species", ab$tx7)
  ab$tx6[which(todo)] <- gsub(" .*$", "", ab$tx7[which(todo)])
  
  # reconstruction of the new full taxonomy from each taxonomic field
  my.ab.result <- data.table::data.table(observation_name=ab$observation_name,  blast_taxonomy= apply(ab, 1, function(r) {
    paste(r[1:7], sep=";", collapse=";")
  }))
  head(my.ab.result)
  dim(my.ab.result)

  # Export of a new taxonomy and a new abundance file
  my.result <- merge(my.mh.result, my.ab.result,by="observation_name",all.y=TRUE)
  final <- merge(my.result, at,by="observation_name",all.y=TRUE)
  final$blast_taxonomy <- apply(final, 1, function(r){
    if (! is.na(r["blast_taxonomy.x"])) {r["blast_taxonomy.x"]} else {r["blast_taxonomy.y"]}
  })
  final$blast_taxonomy.x<-NULL
  final$blast_taxonomy.y<-NULL
  data.table::setcolorder(final, c(2,3,4,5,6,7,8,9,10,1,11:(ncol(final))))
  fordada <- as.data.frame(final)
  fordada <- fordada[order(as.integer(gsub("Cluster_", "", fordada$observation_name))),]
  print("data to find chimeras:")
  print(colnames(fordada))
  # Remove chimeras
  find_chimeras(fordada=fordada, name=name)

  }
#' @examples clean_OTU(abundance_table="Galaxy21-[FROGS_BIOM_to_TSV__abundance.tsv].tsv", affiliation_table="Galaxy22-[FROGS_BIOM_to_TSV__multi-affiliations.tsv].tsv", thresh=0.5, name = "2022JRLAT-16S")