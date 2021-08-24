setwd("/micropollutant_mc2021/ecological group analysis")

library(DESeq2)
library(stringr)
library(gplots)
library(RColorBrewer)
library(data.table)
library(vegan)
library(stats)

mycounts<-read.csv("data for ecological group/even36000.csv")
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
mycounts<-mycounts[,-ncol(mycounts)]
tax<-read.csv('data for ecological group/even36000.csv')
rownames(tax)<-tax[,1]
tax<-tax[,-1]
tax[,'taxonomy']<-as.character(tax[,'taxonomy'])
group<-read.csv("data for ecological group/group.csv")
rownames(group)<-group[,1]
group<-group[,-1]
rrn<-read.csv('data for ecological group/rrn.csv')
rownames(rrn)<-rrn[,1]
rrn<-rrn[,-1]

##control group
counts<-mycounts[,which(group$Treatment=='BEC')]
colData<-group[which(group$Treatment=='BEC'),]
dds <- DESeqDataSetFromMatrix(counts, colData, design= ~ Phase)
dds <- DESeq(dds)
res = results(dds, contrast=c("Phase", "P1", "P2"))
res = res[order(res$pvalue),]
summary(res)
diff_p1p2 <-subset(res,padj < 0.05)
non_diff_p1p2 <-subset(res,padj > 0.05)
res = results(dds, contrast=c("Phase", "P1", "P3"))
res = res[order(res$pvalue),]
diff_p1p3 <-subset(res,padj < 0.05)
non_diff_p1p3 <-subset(res,padj > 0.05)
res = results(dds, contrast=c("Phase", "P2", "P3"))
res = res[order(res$pvalue),]
diff_p2p3 <-subset(res,padj < 0.05)
non_diff_p2p3 <-subset(res,padj > 0.05)

#con.-sensitive
sensitive_otu<-intersect(rownames(subset(diff_p1p2,log2FoldChange > 0)),rownames(subset(diff_p1p3,log2FoldChange > 0)))
sensitive <-  merge(as.data.frame(diff_p1p2[which(diff_p1p2$log2FoldChange > 0&rownames(diff_p1p2)%in%sensitive_otu),]),as.data.frame(diff_p1p3[which(diff_p1p3$log2FoldChange > 0&rownames(diff_p1p3)%in%sensitive_otu),]),by='row.names',sort=FALSE)
rownames(sensitive)<-sensitive$Row.names
sensitive<-sensitive[,c(3,9)]
colnames(sensitive)<-c('log2FoldChange(p1/p2)','log2FoldChange(p1/p3)')
sensitive <-merge(sensitive,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
sensitive[,'category']<-'sensitive'

#con.-opportunistic
opportunistic_otu<-intersect(rownames(subset(diff_p1p2,log2FoldChange < 0)),rownames(subset(diff_p2p3,log2FoldChange > 0)))
opportunistic <-  merge(as.data.frame(diff_p1p2[which(diff_p1p2$log2FoldChange < 0&rownames(diff_p1p2)%in%opportunistic_otu),]),as.data.frame(diff_p2p3[which(diff_p2p3$log2FoldChange > 0&rownames(diff_p2p3)%in%opportunistic_otu),]),by='row.names',sort=FALSE)
rownames(opportunistic)<-opportunistic$Row.names
opportunistic<-opportunistic[,c(3,9)]
colnames(opportunistic)<-c('log2FoldChange(p1/p2)','log2FoldChange(p2/p3)')
opportunistic <-merge(opportunistic,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
opportunistic[,'category']<-'opportunistic'

#con.-tolerant
tolerant_otu<-intersect(rownames(subset(diff_p1p3,log2FoldChange < 0)),rownames(subset(diff_p2p3,log2FoldChange < 0)))
tolerant <-  merge(as.data.frame(diff_p1p3[which(diff_p1p3$log2FoldChange < 0&rownames(diff_p1p3)%in%tolerant_otu),]),as.data.frame(diff_p2p3[which(diff_p2p3$log2FoldChange < 0&rownames(diff_p2p3)%in%tolerant_otu),]),by='row.names',sort=FALSE)
rownames(tolerant)<-tolerant$Row.names
tolerant<-tolerant[,c(3,9)]
colnames(tolerant)<-c('log2FoldChange(p1/p2)','log2FoldChange(p2/p3)')
tolerant <-merge(tolerant,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
tolerant[,'category']<-'tolerant'

sensitive_bec<-sensitive
opportunistic_bec<-opportunistic
tolerant_bec<-tolerant

##BPA
counts<-mycounts[,which(group$Treatment=='BPA')]
colData<-group[which(group$Treatment=='BPA'),]
dds <- DESeqDataSetFromMatrix(counts, colData, design= ~ Phase)
dds <- DESeq(dds)
res = results(dds, contrast=c("Phase", "P1", "P2"))
res = res[order(res$pvalue),]
summary(res)
diff_p1p2 <-subset(res,padj < 0.05)
non_diff_p1p2 <-subset(res,padj > 0.05)
res = results(dds, contrast=c("Phase", "P1", "P3"))
res = res[order(res$pvalue),]
diff_p1p3 <-subset(res,padj < 0.05)
non_diff_p1p3 <-subset(res,padj > 0.05)
res = results(dds, contrast=c("Phase", "P2", "P3"))
res = res[order(res$pvalue),]
diff_p2p3 <-subset(res,padj < 0.05)
non_diff_p2p3 <-subset(res,padj > 0.05)

#BPA-sensitive
sensitive_otu<-intersect(rownames(subset(diff_p1p2,log2FoldChange > 0)),rownames(subset(diff_p1p3,log2FoldChange > 0)))
sensitive <-  merge(as.data.frame(diff_p1p2[which(diff_p1p2$log2FoldChange > 0&rownames(diff_p1p2)%in%sensitive_otu),]),as.data.frame(diff_p1p3[which(diff_p1p3$log2FoldChange > 0&rownames(diff_p1p3)%in%sensitive_otu),]),by='row.names',sort=FALSE)
rownames(sensitive)<-sensitive$Row.names
sensitive<-sensitive[,c(3,7,9,13)]
colnames(sensitive)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p1/p3)','padj2')
sensitive <-merge(sensitive,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
sensitive[,'category']<-'sensitive'

#BPA-opportunistic
opportunistic_otu<-intersect(rownames(subset(diff_p1p2,log2FoldChange < 0)),rownames(subset(diff_p2p3,log2FoldChange > 0)))
opportunistic <-  merge(as.data.frame(diff_p1p2[which(diff_p1p2$log2FoldChange < 0&rownames(diff_p1p2)%in%opportunistic_otu),]),as.data.frame(diff_p2p3[which(diff_p2p3$log2FoldChange > 0&rownames(diff_p2p3)%in%opportunistic_otu),]),by='row.names',sort=FALSE)
rownames(opportunistic)<-opportunistic$Row.names
opportunistic<-opportunistic[,c(3,7,9,13)]
colnames(opportunistic)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p2/p3)','padj2')
opportunistic <-merge(opportunistic,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
opportunistic[,'category']<-'opportunistic'

#BPA-tolerant
tolerant_otu<-intersect(rownames(subset(diff_p1p3,log2FoldChange < 0)),rownames(subset(diff_p2p3,log2FoldChange < 0)))
tolerant <-  merge(as.data.frame(diff_p1p3[which(diff_p1p3$log2FoldChange < 0&rownames(diff_p1p3)%in%tolerant_otu),]),as.data.frame(diff_p2p3[which(diff_p2p3$log2FoldChange < 0&rownames(diff_p2p3)%in%tolerant_otu),]),by='row.names',sort=FALSE)
rownames(tolerant)<-tolerant$Row.names
tolerant<-tolerant[,c(3,7,9,13)]
colnames(tolerant)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p2/p3)','padj2')
tolerant <-merge(tolerant,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
tolerant[,'category']<-'tolerant'

sensitive_BPA<-sensitive
opportunistic_BPA<-opportunistic
tolerant_BPA<-tolerant

##TCS
counts<-mycounts[,which(group$Treatment=='TCS')]
colData<-group[which(group$Treatment=='TCS'),]
dds <- DESeqDataSetFromMatrix(counts, colData, design= ~ Phase)
dds <- DESeq(dds)
res = results(dds, contrast=c("Phase", "P1", "P2"))
res = res[order(res$pvalue),]
summary(res)
diff_p1p2 <-subset(res,padj < 0.05)
non_diff_p1p2 <-subset(res,padj > 0.05)
res = results(dds, contrast=c("Phase", "P1", "P3"))
res = res[order(res$pvalue),]
diff_p1p3 <-subset(res,padj < 0.05)
non_diff_p1p3 <-subset(res,padj > 0.05)
res = results(dds, contrast=c("Phase", "P2", "P3"))
res = res[order(res$pvalue),]
diff_p2p3 <-subset(res,padj < 0.05)
non_diff_p2p3 <-subset(res,padj > 0.05)

#TCS-sensitive
sensitive_otu<-intersect(rownames(subset(diff_p1p2,log2FoldChange > 0)),rownames(subset(diff_p1p3,log2FoldChange > 0)))
sensitive <-  merge(as.data.frame(diff_p1p2[which(diff_p1p2$log2FoldChange > 0&rownames(diff_p1p2)%in%sensitive_otu),]),as.data.frame(diff_p1p3[which(diff_p1p3$log2FoldChange > 0&rownames(diff_p1p3)%in%sensitive_otu),]),by='row.names',sort=FALSE)
rownames(sensitive)<-sensitive$Row.names
sensitive<-sensitive[,c(3,7,9,13)]
colnames(sensitive)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p1/p3)','padj2')
sensitive <-merge(sensitive,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
sensitive[,'category']<-'sensitive'

#TCS-opportunistic
opportunistic_otu<-intersect(rownames(subset(diff_p1p2,log2FoldChange < 0)),rownames(subset(diff_p2p3,log2FoldChange > 0)))
opportunistic <-  merge(as.data.frame(diff_p1p2[which(diff_p1p2$log2FoldChange < 0&rownames(diff_p1p2)%in%opportunistic_otu),]),as.data.frame(diff_p2p3[which(diff_p2p3$log2FoldChange > 0&rownames(diff_p2p3)%in%opportunistic_otu),]),by='row.names',sort=FALSE)
rownames(opportunistic)<-opportunistic$Row.names
opportunistic<-opportunistic[,c(3,7,9,13)]
colnames(opportunistic)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p2/p3)','padj2')
opportunistic <-merge(opportunistic,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
opportunistic[,'category']<-'opportunistic'

#TCS-tolerant
tolerant_otu<-intersect(rownames(subset(diff_p1p3,log2FoldChange < 0)),rownames(subset(diff_p2p3,log2FoldChange < 0)))
tolerant <-  merge(as.data.frame(diff_p1p3[which(diff_p1p3$log2FoldChange < 0&rownames(diff_p1p3)%in%tolerant_otu),]),as.data.frame(diff_p2p3[which(diff_p2p3$log2FoldChange < 0&rownames(diff_p2p3)%in%tolerant_otu),]),by='row.names',sort=FALSE)
rownames(tolerant)<-tolerant$Row.names
tolerant<-tolerant[,c(3,7,9,13)]
colnames(tolerant)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p2/p3)','padj2')
tolerant <-merge(tolerant,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
tolerant[,'category']<-'tolerant'

sensitive_TCS<-sensitive
opportunistic_TCS<-opportunistic
tolerant_TCS<-tolerant

##MI
counts<-mycounts[,which(group$Treatment=='MI')]
colData<-group[which(group$Treatment=='MI'),]
dds <- DESeqDataSetFromMatrix(counts, colData, design= ~ Phase)
dds <- DESeq(dds)
res = results(dds, contrast=c("Phase", "P1", "P2"))
res = res[order(res$pvalue),]
summary(res)
diff_p1p2 <-subset(res,padj < 0.05)
non_diff_p1p2 <-subset(res,padj > 0.05)
res = results(dds, contrast=c("Phase", "P1", "P3"))
res = res[order(res$pvalue),]
diff_p1p3 <-subset(res,padj < 0.05)
non_diff_p1p3 <-subset(res,padj > 0.05)
res = results(dds, contrast=c("Phase", "P2", "P3"))
res = res[order(res$pvalue),]
diff_p2p3 <-subset(res,padj < 0.05)
non_diff_p2p3 <-subset(res,padj > 0.05)

#MI-sensitive
sensitive_otu<-intersect(rownames(subset(diff_p1p2,log2FoldChange > 0)),rownames(subset(diff_p1p3,log2FoldChange > 0)))
sensitive <-  merge(as.data.frame(diff_p1p2[which(diff_p1p2$log2FoldChange > 0&rownames(diff_p1p2)%in%sensitive_otu),]),as.data.frame(diff_p1p3[which(diff_p1p3$log2FoldChange > 0&rownames(diff_p1p3)%in%sensitive_otu),]),by='row.names',sort=FALSE)
rownames(sensitive)<-sensitive$Row.names
sensitive<-sensitive[,c(3,7,9,13)]
colnames(sensitive)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p1/p3)','padj2')
sensitive <-merge(sensitive,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
sensitive[,'category']<-'sensitive'

#MI-opportunistic
opportunistic_otu<-intersect(rownames(subset(diff_p1p2,log2FoldChange < 0)),rownames(subset(diff_p2p3,log2FoldChange > 0)))
opportunistic <-  merge(as.data.frame(diff_p1p2[which(diff_p1p2$log2FoldChange < 0&rownames(diff_p1p2)%in%opportunistic_otu),]),as.data.frame(diff_p2p3[which(diff_p2p3$log2FoldChange > 0&rownames(diff_p2p3)%in%opportunistic_otu),]),by='row.names',sort=FALSE)
rownames(opportunistic)<-opportunistic$Row.names
opportunistic<-opportunistic[,c(3,7,9,13)]
colnames(opportunistic)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p2/p3)','padj2')
opportunistic <-merge(opportunistic,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
opportunistic[,'category']<-'opportunistic'

#MI-tolerant
tolerant_otu<-intersect(rownames(subset(diff_p1p3,log2FoldChange < 0)),rownames(subset(diff_p2p3,log2FoldChange < 0)))
tolerant <-  merge(as.data.frame(diff_p1p3[which(diff_p1p3$log2FoldChange < 0&rownames(diff_p1p3)%in%tolerant_otu),]),as.data.frame(diff_p2p3[which(diff_p2p3$log2FoldChange < 0&rownames(diff_p2p3)%in%tolerant_otu),]),by='row.names',sort=FALSE)
rownames(tolerant)<-tolerant$Row.names
tolerant<-tolerant[,c(3,7,9,13)]
colnames(tolerant)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p2/p3)','padj2')
tolerant <-merge(tolerant,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
tolerant[,'category']<-'tolerant'

sensitive_MI<-sensitive
opportunistic_MI<-opportunistic
tolerant_MI<-tolerant

##MII
counts<-mycounts[,which(group$Treatment=='MII')]
colData<-group[which(group$Treatment=='MII'),]
dds <- DESeqDataSetFromMatrix(counts, colData, design= ~ Phase)
dds <- DESeq(dds)
res = results(dds, contrast=c("Phase", "P1", "P2"))
res = res[order(res$pvalue),]
summary(res)
diff_p1p2 <-subset(res,padj < 0.05)
non_diff_p1p2 <-subset(res,padj > 0.05)
res = results(dds, contrast=c("Phase", "P1", "P3"))
res = res[order(res$pvalue),]
diff_p1p3 <-subset(res,padj < 0.05)
non_diff_p1p3 <-subset(res,padj > 0.05)
res = results(dds, contrast=c("Phase", "P2", "P3"))
res = res[order(res$pvalue),]
diff_p2p3 <-subset(res,padj < 0.05)
non_diff_p2p3 <-subset(res,padj > 0.05)

#MII-sensitive
sensitive_otu<-intersect(rownames(subset(diff_p1p2,log2FoldChange > 0)),rownames(subset(diff_p1p3,log2FoldChange > 0)))
sensitive <-  merge(as.data.frame(diff_p1p2[which(diff_p1p2$log2FoldChange > 0&rownames(diff_p1p2)%in%sensitive_otu),]),as.data.frame(diff_p1p3[which(diff_p1p3$log2FoldChange > 0&rownames(diff_p1p3)%in%sensitive_otu),]),by='row.names',sort=FALSE)
rownames(sensitive)<-sensitive$Row.names
sensitive<-sensitive[,c(3,7,9,13)]
colnames(sensitive)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p1/p3)','padj2')
sensitive <-merge(sensitive,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
sensitive[,'category']<-'sensitive'

#MII-opportunistic
opportunistic_otu<-intersect(rownames(subset(diff_p1p2,log2FoldChange < 0)),rownames(subset(diff_p2p3,log2FoldChange > 0)))
opportunistic <-  merge(as.data.frame(diff_p1p2[which(diff_p1p2$log2FoldChange < 0&rownames(diff_p1p2)%in%opportunistic_otu),]),as.data.frame(diff_p2p3[which(diff_p2p3$log2FoldChange > 0&rownames(diff_p2p3)%in%opportunistic_otu),]),by='row.names',sort=FALSE)
rownames(opportunistic)<-opportunistic$Row.names
opportunistic<-opportunistic[,c(3,7,9,13)]
colnames(opportunistic)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p2/p3)','padj2')
opportunistic <-merge(opportunistic,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
opportunistic[,'category']<-'opportunistic'

#MII-tolerant
tolerant_otu<-intersect(rownames(subset(diff_p1p3,log2FoldChange < 0)),rownames(subset(diff_p2p3,log2FoldChange < 0)))
tolerant <-  merge(as.data.frame(diff_p1p3[which(diff_p1p3$log2FoldChange < 0&rownames(diff_p1p3)%in%tolerant_otu),]),as.data.frame(diff_p2p3[which(diff_p2p3$log2FoldChange < 0&rownames(diff_p2p3)%in%tolerant_otu),]),by='row.names',sort=FALSE)
rownames(tolerant)<-tolerant$Row.names
tolerant<-tolerant[,c(3,7,9,13)]
colnames(tolerant)<-c('log2FoldChange(p1/p2)','padj1','log2FoldChange(p2/p3)','padj2')
tolerant <-merge(tolerant,as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
tolerant[,'category']<-'tolerant'

sensitive_MII<-sensitive
opportunistic_MII<-opportunistic
tolerant_MII<-tolerant

##combine sensitive, opportunistic, tolerant OTUs into one dataframe and assign with taxonomy
sensitive_all<-rbind(data.frame(Treatment='BPA',sensitive_BPA[,c(1:5,27)]),data.frame(Treatment='TCS',sensitive_TCS[,c(1:5,27)]),data.frame(Treatment='MI',sensitive_MI[,c(1:5,27)]),data.frame(Treatment='MII',sensitive_MII[,c(1:5,27)]))
opportunistic_all<-rbind(data.frame(Treatment='BPA',opportunistic_BPA[,c(1:5,27)]),data.frame(Treatment='TCS',opportunistic_TCS[,c(1:5,27)]),data.frame(Treatment='MI',opportunistic_MI[,c(1:5,27)]))
tolerant_all<-rbind(data.frame(Treatment='BPA',tolerant_BPA[,c(1:5,27)]),data.frame(Treatment='TCS',tolerant_TCS[,c(1:5,27)]),data.frame(Treatment='MI',tolerant_MI[,c(1:5,27)]),data.frame(Treatment='MII',tolerant_MII[,c(1:5,27)]))

colnames(sensitive_all)[2]<-'OTUID'
colnames(opportunistic_all)[2]<-'OTUID'
colnames(tolerant_all)[2]<-'OTUID'

for(i in 1:nrow(sensitive_all)){
  sensitive_all[i,'taxonomy']<-tax[which(rownames(tax)==sensitive_all[i,'OTUID']),'taxonomy']
}

for(i in 1:nrow(opportunistic_all)){
  opportunistic_all[i,'taxonomy']<-tax[which(rownames(tax)==opportunistic_all[i,'OTUID']),'taxonomy']
}

for(i in 1:nrow(tolerant_all)){
  tolerant_all[i,'taxonomy']<-tax[which(rownames(tax)==tolerant_all[i,'OTUID']),'taxonomy']
}

#export tables
write.csv(sensitive_all,'sensitive_all.csv')
write.csv(opportunistic_all,'opportunistic_all.csv')
write.csv(tolerant_all,'tolerant_all.csv')

#calculate relative abundance of OTU per batch
all_otu<-rbind(sensitive_all[,c(1:2,7)],opportunistic_all[,c(1:2,7)],tolerant_all[,c(1:2,7)])
all_otu<-data.frame(all_otu,P1=rep(NA,nrow(all_otu)),P2=rep(NA,nrow(all_otu)),P3=rep(NA,nrow(all_otu)))
for (i in 1:nrow(all_otu)){
  ra_otu<-mycounts[which(rownames(mycounts)==all_otu[i,'OTUID']),]
  ra_otu<-ra_otu[grep(all_otu[i,'Treatment'],colnames(mycounts))]
  ra_otu1<-cbind(ra_otu[grep('I1',colnames(ra_otu))],ra_otu[grep('I2',colnames(ra_otu))])
  all_otu[i,'P1']<-mean(as.numeric(ra_otu1))/36000
  ra_otu2<-cbind(ra_otu[grep('I3',colnames(ra_otu))],ra_otu[grep('I4',colnames(ra_otu))])
  all_otu[i,'P2']<-mean(as.numeric(ra_otu2))/36000
  ra_otu3<-cbind(ra_otu[grep('I5',colnames(ra_otu))],ra_otu[grep('I6',colnames(ra_otu))],ra_otu[grep('I7',colnames(ra_otu))])
  all_otu[i,'P3']<-mean(as.numeric(ra_otu3))/36000
  all_otu$I1[i]<-mean(as.numeric(ra_otu[grep('I1',colnames(ra_otu))]))/36000
  all_otu$I2[i]<-mean(as.numeric(ra_otu[grep('I2',colnames(ra_otu))]))/36000
  all_otu$I3[i]<-mean(as.numeric(ra_otu[grep('I3',colnames(ra_otu))]))/36000
  all_otu$I4[i]<-mean(as.numeric(ra_otu[grep('I4',colnames(ra_otu))]))/36000
  all_otu$I5[i]<-mean(as.numeric(ra_otu[grep('I5',colnames(ra_otu))]))/36000
  all_otu$I6[i]<-mean(as.numeric(ra_otu[grep('I6',colnames(ra_otu))]))/36000
  all_otu$I7[i]<-mean(as.numeric(ra_otu[grep('I7',colnames(ra_otu))]))/36000
}

#with relative abundance per phase
all_otu2<-melt(all_otu[1:6],id=colnames(all_otu)[1:3])
all_otu3<-all_otu2
all_otu3<-all_otu2[which(all_otu2$category!='opportunistic'),]#'opportunistic' group excluded
#with relative abundance per batch
all_otu_batch<-melt(all_otu[c(1:3,7:13)],id=colnames(all_otu)[1:3])
all_otu_b<-all_otu_batch
all_otu_b<-all_otu_b[which(all_otu_b$category!='opportunistic'),]#'opportunistic' group excluded

#Fig6a heatmap
my_palette<-colorRampPalette(c("royalblue4","blue","white","pink","tomato"))((n=80))
breaks <- seq(from = 0, to = max(all_otu3$value), length = 81)

#bpa
bpa_otu_b<-all_otu_b[which(all_otu_b$Treatment=='BPA'),]
bpa_otu_b$OTUID<-factor(bpa_otu_b$OTUID,levels = bpa_otu_b[which(bpa_otu_b$variable=='I1'),'OTUID'])
bpa_otu_b<-dcast(bpa_otu_b, Treatment+OTUID+category~variable)
rownames(bpa_otu_b)<-bpa_otu_b$OTUID

par(family = 'Times') 
fig6a_bpa_b<-heatmap.2(as.matrix(bpa_otu_b[,c(4:10)]),col=my_palette, dendrogram="none",Rowv=NA,Colv=NA, density.info="none",trace="none", 
                   key=TRUE,scale='row',keysize=1, srtCol = 0, cexCol=1.2, cexRow=0.8,adjRow=c(0,0), adjCol=c(0.5,0.5),lmat=rbind(c(5,0,4), c(3,2,1)), lhei=c(1,6),lwid=c(2,3,0.3),
                   margins =c(2,4),
                   RowSideColors = c(rep("#006600", length(which(bpa_otu_b$category=='sensitive'))), rep("orange", length(which(bpa_otu_b$category=='opportunistic'))),rep("#00CCFF", length(which(bpa_otu_b$category=='tolerant')))))
dev.off()

#tcs
tcs_otu_b<-all_otu_b[which(all_otu_b$Treatment=='TCS'),]
tcs_otu_b$OTUID<-factor(tcs_otu_b$OTUID,levels = tcs_otu_b[which(tcs_otu_b$variable=='I1'),'OTUID'])
tcs_otu_b<-dcast(tcs_otu_b, Treatment+OTUID+category~variable)
rownames(tcs_otu_b)<-tcs_otu_b$OTUID

par(family = 'Times') 
fig6a_tcs_b<-heatmap.2(as.matrix(tcs_otu_b[,c(4:10)]),col=my_palette, dendrogram="none",Rowv=NA,Colv=NA, density.info="none",trace="none", 
                   key=TRUE,scale='row',keysize=1, srtCol = 0, cexCol=1.2, cexRow=0.8,adjRow=c(0,0), adjCol=c(0.5,0.5),lmat=rbind(c(5,0,4), c(3,2,1)), lhei=c(1,6),lwid=c(2,3,0.3),
                   margins =c(2,4),
                   RowSideColors = c(rep("#006600", length(which(tcs_otu_b$category=='sensitive'))), rep("orange", length(which(tcs_otu_b$category=='opportunistic'))),rep("#00CCFF", length(which(tcs_otu_b$category=='tolerant')))))
dev.off()

#mi
mi_otu_b<-all_otu_b[which(all_otu_b$Treatment=='MI'),]
mi_otu_b$OTUID<-factor(mi_otu_b$OTUID,levels = mi_otu_b[which(mi_otu_b$variable=='I1'),'OTUID'])
mi_otu_b<-dcast(mi_otu_b, Treatment+OTUID+category~variable)
rownames(mi_otu_b)<-mi_otu_b$OTUID

par(family = 'Times') 
p_mi_b<-heatmap.2(as.matrix(mi_otu_b[,c(4:10)]),col=my_palette, dendrogram="none",Rowv=NA,Colv=NA, density.info="none",trace="none", 
                  key=TRUE,scale='row',keysize=1, srtCol = 0, cexCol=1.2, cexRow=0.8,adjRow=c(0,0), adjCol=c(0.5,0.5),lmat=rbind(c(5,0,4), c(3,2,1)), lhei=c(1,6),lwid=c(2,3,0.3),
                  margins =c(2,4),
                  RowSideColors = c(rep("#006600", length(which(mi_otu_b$category=='sensitive'))), rep("orange", length(which(mi_otu_b$category=='opportunistic'))),rep("#00CCFF", length(which(mi_otu_b$category=='tolerant')))))
dev.off()

#mii
mii_otu_b<-all_otu_b[which(all_otu_b$Treatment=='MII'),]
mii_otu_b$OTUID<-factor(mii_otu_b$OTUID,levels = mii_otu_b[which(mii_otu_b$variable=='I1'),'OTUID'])
mii_otu_b<-dcast(mii_otu_b, Treatment+OTUID+category~variable)
rownames(mii_otu_b)<-mii_otu_b$OTUID

par(family = 'Times') 
p_mii_b<-heatmap.2(as.matrix(mii_otu_b[,c(4:10)]),col=my_palette, dendrogram="none",Rowv=NA,Colv=NA, density.info="none",trace="none", 
                   key=TRUE,scale='row',keysize=1, srtCol = 0, cexCol=1.2, cexRow=0.8,adjRow=c(0,0), adjCol=c(0.5,0.5),lmat=rbind(c(5,0,4), c(3,2,1)), lhei=c(1,6),lwid=c(2,3,0.3),
                   margins =c(2,4),
                   RowSideColors = c(rep("#006600", length(which(mii_otu_b$category=='sensitive'))), rep("orange", length(which(mii_otu_b$category=='opportunistic'))),rep("#00CCFF", length(which(mii_otu_b$category=='tolerant')))))
dev.off()

#cluster
#by rrn
#without control
hc_rrn_d=hclust(vegdist(t(rrn[4:15,]),method="bray"), method = "complete")

plot(hc_rrn_d,labels =c('B1','B2','B3','B4','B5','B6','B7'), main = "Cluster Dendrogram",ylab =NA,xlab = NA,sub= NA,family = "Times",cex.axis=1)

