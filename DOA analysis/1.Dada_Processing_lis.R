#this manuscript is modified from https://doi.org/10.5281/zenodo.3817698. 
#Vila, J. C., Liu, Y. Y., & Sanchez, A. (2020). Dissimilarity-Overlap analysis of replicate enrichment communities. The ISME Journal, 14(10), 2505-2513.

library(data.table)
library(readr)
set.seed(1)
#Script 1 for Extract data  from DADA2 output and rarefies to constant read depth. 
# Stores into a melted data.frame with standardized columns for subsequent analysis.

rarefy <- function(dat,n =min(colSums(dat))){
  # normalize to sample count
  for(i in 1:ncol(dat)){dat
    old_column = dat[,i,with=FALSE][[1]]
    elements = factor(row.names(dat),levels=row.names(dat))
    sample = table(sample(rep(elements,old_column),n,replace =FALSE))
    dat[,(i) := as.numeric(sample)]
  }
  return(dat)
}

assign_taxonomy_id <- function(dat){
  #Asigns ID to the higest available taxonomic level
  dat$ESV_ID = dat$Genus
  dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Family[which(is.na(dat$ESV_ID))]
  dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Order[which(is.na(dat$ESV_ID))]
  dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Class[which(is.na(dat$ESV_ID))]
  dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Phylum[which(is.na(dat$ESV_ID))]
  dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Kingdom[which(is.na(dat$ESV_ID))]
  dat$ESV_ID = as.factor(make.unique(dat$ESV_ID))
  #Make NA use ESV_ID to avoid counting them as the same group
  dat[is.na(dat$Genus)]$Genus = as.character(dat[is.na(dat$Genus)]$ESV_ID)
  dat[is.na(dat$Family)]$Family = dat[is.na(dat$Family)]$Genus
  dat[is.na(dat$Order)]$Order = dat[is.na(dat$Order)]$Family
  dat[is.na(dat$Class)]$Class = dat[is.na(dat$Class)]$Order
  dat[is.na(dat$Phylum)]$Phylum = dat[is.na(dat$Phylum)]$Class
  
  return(dat)
}
set.seed(1) # To ensure reproducibility

#Load Data
setwd("/micropollutant_mc2021/DOC analysis/data for DOC")
Taxonomy_Data = fread('taxonomy.csv')
Abundance_Data = fread('otu_table.csv') #Data is actually at an ESV level
setwd("/micropollutant_mc2021/DOC analysis")

#Rarefy
Abundance_Data = as.matrix(Abundance_Data[,2:ncol(Abundance_Data)])
Abundance_Data = data.table(Abundance_Data)
Abundance_Data = rarefy(Abundance_Data) 

Taxonomy_Data = assign_taxonomy_id(Taxonomy_Data) #Assign Taxonomy ID
Abundance_Data = t(t(Abundance_Data)/colSums(Abundance_Data)) #Calculated Relative Abundance of ESVs

Abundance_Data = data.table(Abundance_Data)
Abundance_Data$ESV_ID = Taxonomy_Data$ESV_ID

#Now convert matrix into a data.frame using melt
Abundance_Data2 = melt(Abundance_Data,id= c('ESV_ID'),variable.name ='Sample_ID',value.name = 'Relative_Abundance')
Abundance_Data<-Abundance_Data2
Abundance_Data = Abundance_Data[Relative_Abundance>0,]

#Extract batch, treatment,replicate number from sample id
Abundance_Data$batch =sapply(Abundance_Data$Sample_ID,function(x) substr(as.character(x),1,2))
Abundance_Data$treatment =sapply(Abundance_Data$Sample_ID,function(x) substr(as.character(x),3,5))
Abundance_Data$replicate =sapply(Abundance_Data$Sample_ID,function(x) substr(as.character(x),8,8))
Abundance_Data$ESV_ID = droplevels(Abundance_Data$ESV_ID)
#correction
Abundance_Data[which(Abundance_Data$batch=='I0'&Abundance_Data$treatment=='1MI'),'batch']<-'I1'
Abundance_Data[which(Abundance_Data$batch=='I1'&Abundance_Data$treatment=='1MI'),'treatment']<-'MI'
Abundance_Data[which(Abundance_Data$batch=='I0'&Abundance_Data$treatment=='2MI'),'batch']<-'I2'
Abundance_Data[which(Abundance_Data$batch=='I2'&Abundance_Data$treatment=='2MI'),'treatment']<-'MI'
Abundance_Data[which(Abundance_Data$batch=='I0'&Abundance_Data$treatment=='3MI'),'batch']<-'I3'
Abundance_Data[which(Abundance_Data$batch=='I3'&Abundance_Data$treatment=='3MI'),'treatment']<-'MI'
Abundance_Data[which(Abundance_Data$batch=='I0'&Abundance_Data$treatment=='4MI'),'batch']<-'I4'
Abundance_Data[which(Abundance_Data$batch=='I4'&Abundance_Data$treatment=='4MI'),'treatment']<-'MI'
Abundance_Data[which(Abundance_Data$batch=='I0'&Abundance_Data$treatment=='5MI'),'batch']<-'I5'
Abundance_Data[which(Abundance_Data$batch=='I5'&Abundance_Data$treatment=='5MI'),'treatment']<-'MI'
Abundance_Data[which(Abundance_Data$batch=='I0'&Abundance_Data$treatment=='6MI'),'batch']<-'I6'
Abundance_Data[which(Abundance_Data$batch=='I6'&Abundance_Data$treatment=='6MI'),'treatment']<-'MI'
Abundance_Data[which(Abundance_Data$batch=='I0'&Abundance_Data$treatment=='7MI'),'batch']<-'I7'
Abundance_Data[which(Abundance_Data$batch=='I7'&Abundance_Data$treatment=='7MI'),'treatment']<-'MI'

#Merge with taxonomy
merged_data = merge(Taxonomy_Data,Abundance_Data)

#Save Eveything
fwrite(merged_data,'Emergent_Comunity_Data.csv')

