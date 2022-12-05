# date : 20211208
# version : 0.01
# authors : Sylvain Piry (CBGP)
# licence : GPL

rm(list=ls())# Reset

#####################
##                 ##
##      SHORT      ##
##   DESCRIPTION   ##
##                 ##
#####################

## Work on raw abundance and multihits text files from Frogs.

## Correct the taxonomy of each cluster of the abundance tsv file tagged 'multi-affiliation', from decision rules applied on the cleaned multihits tsv file:

# Clean the taxonomy (7 ranks: Kingdom,Phylum,Class,Order,Family,Genus,Species) of each hit of the two files for:
# - the special characters '[]()_"0123456789.' in all 7 taxonomic ranks.
# - 'mitochondria' (in the Family rank) and 'chloroplast' (in the Class rank), which are replaced by plantae in all 7 taxonomic ranks.
# - 'unidentified', 'sp.' and 'bacterium', which are replaced by 'unknown species' in the Species rank .
# - unexpected number of words in all 7 taxonomic ranks.
# - for the species and genus ranks only: if the number of 'unknown species' hits < % threshold then ignore them.
# - for all taxonomic ranks: 
  # - if there is a unique taxonomic name, then rename by this taxonomic name.
  # - if there is a unique taxonomic name and an 'unknown', then rename by 'unknown'.
  # - if there are two or more taxonomic names, then name 'multi-affiliation'.

## Identify chimeric sequences using the isBimeraDenovo(...) function in the dada2 library

## Export an abundance table tsv file cleaned for taxonomy and filtered for chimeras.

###################
##               ##
##    INPUT      ##
## REQUIREMENTS  ##
##               ##
###################

## In the set folder: a multihit tsv file and an abundance tsv file.
## Search for "to modify" in this script and change values accordingly to your data.
## Dependencies: packages to install are listed below.
library(data.table)
library(tidyr)
library(dplyr)
library(dada2) # BioConductor

###################
##               ##
##  DATA IMPORT  ##
##               ##
###################

## Load input files from data folder
setwd("/home/benoitla/Documents/_Analyses/2022_Frago/") #to modify
name<-"Frago2022-16S" #name will be used to import files and write names before extension of the different output tables. to modify

mt<-fread(paste(name,"_multihits.txt",sep=""),sep="\t",header=TRUE) # to modify (name of the multihits tsv file)
mh <- data.table(bidouille=mt$blast_taxonomy, blast_taxonomy=mt$blast_taxonomy, observation_name=mt$`#observation_name`)
rm(mt)

at<-fread(paste(name,"_abundance.txt",sep=""),sep="\t",header=TRUE)# to modify (name of the abundance file)
ab <- data.table(bidouille=at$`blast_taxonomy`, blast_taxonomy=at$`blast_taxonomy`, observation_name=at$observation_name)

###################
##               ##
##  CLEANING OF  ##
##      THE      ##
##   MULTIHITS   ##
##     FILE      ##
##               ##
###################

## Cleaning for the special characters []()_"0123456789.
mh$bidouille<-gsub("[\\[\\]()_\"0-9\\.]","",mh$bidouille,perl=TRUE)
mh$bidouille<-gsub(" +;",";",mh$bidouille,perl=TRUE)

## Deconstruction of the taxonomic column ('blast_taxonomy') in 7 taxonomic fields
toto<- mh[ , "bidouille"] %>% separate(bidouille,c("tx1","tx2","tx3","tx4","tx5","tx6","tx7"),";")
mh<-cbind(toto,mh[, c("observation_name" ,"blast_taxonomy")]) # observation name is the cluster ID and useful to link with the abundance tsv file
rm(toto)

## Replacement of 'chloroplast' and 'mitochondria' by 'plantae' in all 7 taxonomic fields
mh[grepl("Chloroplast",mh$tx4,ignore.case=TRUE)|grepl("Chloroplast",mh$tx4,ignore.case=TRUE)|grepl("Mitochondria",mh$tx5,ignore.case=TRUE),1:7]<-c("Plantae")
#rep("Plantae",7) ne marche plus

## Overwriting of the species field with 'unknown species' if presence of 'incertae sedis' in any of upper taxonomic fields 
## Overwriting of the genus field with 'unknown genus'  if presence of 'incertae sedis' in any of upper taxonomic fields
mh[!(grepl("incertae sedis",mh$tx2,ignore.case=TRUE)|grepl("incertae sedis",mh$tx3,ignore.case=TRUE)|grepl("incertae sedis",mh$tx4,ignore.case=TRUE)|grepl("incertae sedis",mh$tx5,ignore.case=TRUE)|grepl("incertae sedis",mh$tx6,ignore.case=TRUE))&grepl("incertae sedis",mh$tx7,ignore.case=TRUE),7]<-c("unknown species")
mh[grepl("incertae sedis",mh$tx2,ignore.case=TRUE)|grepl("incertae sedis",mh$tx3,ignore.case=TRUE)|grepl("incertae sedis",mh$tx4,ignore.case=TRUE)|grepl("incertae sedis",mh$tx5,ignore.case=TRUE)|grepl("incertae sedis",mh$tx6,ignore.case=TRUE),6:7]<-list("unknown genus","unknown species")

## Cleaning for more than 1 word in the 6 first taxonomic fields (except if 'unknown' is found)
### Prior: removal of 'Candidatus' in each field
mh <- as.data.frame(mh)
for(j in 1:6) {
  mh[,j] <- gsub("^Candidatus ", "", mh[,j])
  my.stuff <- as.data.frame(mh[,j])
  tf <- !grepl("unknown",my.stuff[,1],ignore.case=TRUE)
  mh[which(tf),j] <- gsub(" +.*$","", mh[which(tf),j], perl=TRUE)
}
mh <- as.data.table(mh)

## Cleaning for more than 2 words in the 7th taxonomic field
## Prior: removal of 'Candidatus' in the species field
mh$tx7<-gsub("^Candidatus ","",mh$tx7,perl=TRUE)
mh$tx7<-gsub(" +$", "", gsub("^(([^ ]+ ?){1,2}).*$","\\1",mh$tx7,perl=TRUE),perl=TRUE)

## Replacement in the species field of 'unidentified' or 'metagenome' anywhere in the field by 'unknown species' 
## Replacement in the species field of 'sp' at the end of the field by 'unknown species'
## Replacement in the species field of 'bacterium' or 'uncultured' at the beginning of the field by 'unknown species'
## Replacement of the species field by 'unknown species' if empty
mh$tx7<-gsub("^.*(unidentified|metagenome).*$","unknown species",mh$tx7,perl=TRUE)
mh$tx7<-gsub("^.*(sp|eae bacterium|ales bacterium)$","unknown species",mh$tx7,perl=TRUE)
mh$tx7<-gsub("^(uncultured|bacterium).*$","unknown species",mh$tx7,perl=TRUE)
mh$tx7[grepl("^ *$",mh$tx7)]<-"unknown species"

## Overwriting of the species field with 'unknown species' if it begins by a lower case 
todo <- grepl("^[a-z]", mh$tx7) & !grepl("unknown species", mh$tx7)
mh$tx7[which(todo)] <- "unknown species"

## Overwriting of the genus field  (if unknown) with the first word of the species field (if not unknown) if empty 
todo <- mh$tx6=="unknown genus" & grepl(" ", mh$tx7) & !grepl("unknown species", mh$tx7)
mh$tx6[which(todo)] <- gsub(" .*$", "", mh$tx7[which(todo)])

## Ignoring of 'unknown species' or 'unknown genus' under a threshold
thresh=0.5 # to modify % threshold below which 'unknown' are ignored
mh$clus <- as.integer(gsub("Cluster_","",mh$observation_name))
mh$blast_taxonomy<-NULL

my.genus.result <- as.data.table(unique(mh %>% add_count(clus, tx1, tx2, tx3, tx4, tx5, tx6, observation_name)))
my.genus.result$ntot<-as.data.table(unique(mh %>% add_count(observation_name)))$n
setkey(my.genus.result, clus, n)
my.genus.result$prop.genus <-my.genus.result$n/my.genus.result$ntot
my.genus.result$ignore.genus <-FALSE
my.genus.result$ignore.genus <- my.genus.result$ignore.genus | (grepl("unknown",my.genus.result$tx6,perl=TRUE,ignore.case=TRUE) & my.genus.result$prop.genus <thresh)
#nrow(my.genus.result)
#names(my.genus.result)

my.species.result <- as.data.table(unique(mh %>% add_count(clus, tx1, tx2, tx3, tx4, tx5, tx6, tx7, observation_name)))
my.species.result$ntot<-as.data.table(unique(mh %>% add_count(observation_name)))$n
setkey(my.species.result, clus, n)
my.species.result$prop.species<-my.species.result$n/my.species.result$ntot
my.species.result$ignore.species<-FALSE
my.species.result$ignore.species <- my.species.result$ignore.species | (grepl("unknown",my.species.result$tx7,perl=TRUE,ignore.case=TRUE) & my.species.result$prop.species<thresh)
#nrow(my.species.result)
#names(my.species.result)

my.first.result<-merge(my.genus.result, my.species.result, by = names(my.species.result)[1:8])[,c(1:8,13,18)]

my.second.result <- data.table(observation_name=unique(my.first.result$observation_name))
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
my.mh.result <- data.table(observation_name=my.second.result$observation_name,  blast_taxonomy= apply(my.second.result, 1, function(r) {
  paste(r[2:8], sep=";", collapse=";")
}))
#nrow(my.mh.result)

###################
##               ##
##  CLEANING OF  ##
##      THE      ##
##   ABUNDANCE   ##
##     FILE      ##
##               ##
###################

## Cleaning for some special characters
ab$bidouille<-gsub("[\\[\\]()_\"0-9\\.]","",ab$bidouille,perl=TRUE)
ab$bidouille<-gsub(" +;",";",ab$bidouille,perl=TRUE)

## Deconstruction of the taxonomic column in 7 taxonomic fields
toto<- ab[ , "bidouille"] %>% separate(bidouille,c("tx1","tx2","tx3","tx4","tx5","tx6","tx7"),";")
ab<-cbind(toto,ab[, c("observation_name" ,"blast_taxonomy")])
rm(toto)

## Replacement of 'chloroplast' and 'mitochondria' by 'plantae' in all 7 taxonomic fields
ab[grepl("Chloroplast",ab$tx4,ignore.case=TRUE)|grepl("Mitochondria",ab$tx5,ignore.case=TRUE),1:7]<-c("Plantae")

## Overwriting of the species field with 'unknown species' if presence of 'incertae sedis' in any of upper taxonomic fields 
## Overwriting of the genus field with 'unknown genus'  if presence of 'incertae sedis' in any of upper taxonomic fields
ab[!(grepl("incertae sedis",ab$tx2,ignore.case=TRUE)|grepl("incertae sedis",ab$tx3,ignore.case=TRUE)|grepl("incertae sedis",ab$tx4,ignore.case=TRUE)|grepl("incertae sedis",ab$tx5,ignore.case=TRUE)|grepl("incertae sedis",ab$tx6,ignore.case=TRUE))&grepl("incertae sedis",ab$tx7,ignore.case=TRUE),7]<-c("unknown species")
ab[grepl("incertae sedis",ab$tx2,ignore.case=TRUE)|grepl("incertae sedis",ab$tx3,ignore.case=TRUE)|grepl("incertae sedis",ab$tx4,ignore.case=TRUE)|grepl("incertae sedis",ab$tx5,ignore.case=TRUE)|grepl("incertae sedis",ab$tx6,ignore.case=TRUE),6:7]<-list("unknown genus","unknown species")

## Cleaning for more than 1 word in the 6 first taxonomic fields (except if 'unknown' is found)
### Prior: removal of 'Candidatus' in each field
ab <- as.data.frame(ab)
for(j in 1:6) {
  ab[,j] <- gsub("^Candidatus ", "", ab[,j])
  my.stuff <- as.data.frame(ab[,j])
  tf <- !grepl("unknown",my.stuff[,1],ignore.case=TRUE)
  ab[which(tf),j] <- gsub(" +.*$","", ab[which(tf),j], perl=TRUE)
}
ab <- as.data.table(ab)

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
my.ab.result <- data.table(observation_name=ab$observation_name,  blast_taxonomy= apply(ab, 1, function(r) {
  paste(r[1:7], sep=";", collapse=";")
}))

########################
##                    ##
##  WRITING OF A NEW  ##
##  TAXONOMY IN A NEW ##
##   ABUNDANCE FILE   ##
##                    ##
########################

my.result <- merge(my.mh.result, my.ab.result,by="observation_name",all.y=TRUE)
final <- merge(my.result, at,by="observation_name",all.y=TRUE)
final$blast_taxonomy <- apply(final, 1, function(r){
  if (! is.na(r["blast_taxonomy.x"])) {r["blast_taxonomy.x"]} else {r["blast_taxonomy.y"]}
})
final$blast_taxonomy.x<-NULL
final$blast_taxonomy.y<-NULL
setcolorder(final,c(2,3,4,5,6,7,8,9,10,1,11:(ncol(final))))
fordada<-final[order(as.integer(gsub("Cluster_","",final$observation_name)))]
fordada<-as.data.frame(fordada)

####################
##                ##
##     REMOVE     ##
##    CHIMERAS    ##
##                ##    
####################

## Tune data structure
samples <- colnames(fordada)[seq(1+which(colnames(fordada)=="observation_sum"), ncol(fordada))]
fordada.m <- data.matrix(fordada[ , seq(1+which(colnames(fordada)=="observation_sum"), ncol(fordada))])
rownames(fordada.m) <- fordada$seed_sequence
fordada.t <- t(fordada.m)

##
bimeras.v <- isBimeraDenovo(getUniques(fordada.t), minFoldParentOverAbundance=4, verbose=TRUE)
bimeras.df <- data.frame(x=as.logical(bimeras.v), seq=names(bimeras.v))
write.table(fordada[bimeras.df$x==FALSE,], file = paste(name,"-cleaned_abundance.txt",sep=""), sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)#writing of filtered abundance table (without identified chimeras)
write.table(fordada[bimeras.df$x==TRUE,], file = paste(name,"_chimeras_list.txt",sep=""), sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)#downloading of identified chimeras

## Summary file
number_of_total_clusters = nrow(fordada)
number_of_bimeras = nrow(fordada[bimeras.df$x==FALSE,])
number_of_nobimeras = nrow(fordada[bimeras.df$x==TRUE,])
summary.df = data.frame(number_of_total_clusters, number_of_bimeras, number_of_nobimeras)
write.table(summary.df,file=paste(name,"_chimeras_summary.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

