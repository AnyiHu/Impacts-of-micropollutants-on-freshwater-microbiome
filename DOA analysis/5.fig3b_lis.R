#this manuscript is modified from https://doi.org/10.5281/zenodo.3817698. 
#Vila, J. C., Liu, Y. Y., & Sanchez, A. (2020). Dissimilarity-Overlap analysis of replicate enrichment communities. The ISME Journal, 14(10), 2505-2513.

setwd("/micropollutant_mc2021/DOC analysis")
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
windowsFonts(Times=windowsFont("Times New Roman"))
#rm(list=ls())

lab_tr<-c('BEC'='Con.','BPA'='BPA','TCS'='TCS','MI'='MI','MII'='MII')

plot_lowess<- function(dat,fn,med){
  #Plot DOC using lowess curve fitting approach
  #med is the median overlap used,
  #fn is the file contained the boostrapped curve dta
  sdata= fread(fn)
  m = loess(dat$Dissimilarity ~ dat$Overlap,span=0.2,family="symmetric",iterations=5,statistics='none')
  sdata$y = predict(m,sdata$xs)
  sdata = sdata[!is.na(sdata$y)]
  p<- ggplot(dat) +
    geom_point(data=dat,aes(x= Overlap,y = Dissimilarity),size=1,alpha=0.05) +
    geom_segment(x=med,xend=med,y=0,yend=sqrt(log(2)),linetype=2,col='Red') +
    geom_ribbon(data=sdata,aes(x=xs,ymin=LCI_lowess,ymax=UCI_lowess),alpha=0.5,fill='Red') + 
    geom_line(data=sdata,aes(x=xs,y=y),col='Red',size=2) +
    scale_x_continuous(breaks=c(0,1),limits=c(0,1)) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.83))  + theme(axis.line = element_line(size=1)) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title=element_text(size=16))
  return(p)
}


plot_olm<- function(dat,fn,med){
  #Plot olm over subset of data with above med overlap
  # Fn is file containing bootstrapped data
  m = lm(Dissimilarity ~ Overlap, dat[Overlap>med])
  slope =floor(signif(m$coefficients[2],2)*100)/100
  sdata= fread(fn)
  dat = dat[Overlap>med,]
  sdata = sdata[xs>min(dat$Overlap),]
  pvalue = mean(sdata$PValue)
  if(pvalue ==0){pvalue = 'p < 0.002'}else{
    pvalue = paste('p = ',signif(pvalue,2),sep='')
  }
  
  p<- ggplot(dat) +
    geom_point(data=dat,aes(x= Overlap,y = Dissimilarity),size=1,alpha=0.05) +
    theme_pubr() + 
    geom_line(data=sdata,aes(x=xs,y=mean_lm),col='Blue',size=2) +
    geom_ribbon(data=sdata,aes(x=xs,ymin=LCI_lm,ymax=UCI_lm),alpha=0.5,fill='Blue') +
    scale_x_continuous(breaks=c(floor(med*100)/100,1),limits=c(floor(med*100)/100,1)) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9)) + theme(axis.line = element_line(size=1)) +
    annotate(geom='text',x = floor(min(dat$Overlap)*100)/100,y=0.9,label=paste('\t\t\t\t m =', slope,' ' , pvalue),size=2) +
    theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  return(p)
}

plot_olm_cs<- function(dat,med){
  #Plot olm over subset of data with above med overlap gationg by carbon source for figure C and supplementary
  # Fn is file containing bootstrapped data
  #Subset
  bec = dat[treatment_1 == 'BEC']
  bpa = dat[treatment_1 == 'BPA']
  tcs = dat[treatment_1 == 'TCS']
  mi = dat[treatment_1 == 'MI']
  mii = dat[treatment_1 == 'MII']
  #media
  mbec = lm(Dissimilarity ~ Overlap, bec[Overlap>med])
  mbpa =lm(Dissimilarity ~ Overlap, bpa[Overlap>med])
  mtcs =lm(Dissimilarity ~ Overlap, tcs[Overlap>med])
  mmi =lm(Dissimilarity ~ Overlap, mi[Overlap>med])
  mmii =lm(Dissimilarity ~ Overlap, mii[Overlap>med])
  
  slopebec =floor(signif(mbec$coefficients[2],2)*100)/100
  slopebpa =floor(signif(mbpa$coefficients[2],2)*100)/100
  slopetcs =floor(signif(mtcs$coefficients[2],2)*100)/100
  slopemi =floor(signif(mmi$coefficients[2],2)*100)/100
  slopemii =floor(signif(mmii$coefficients[2],2)*100)/100
  
  sdatabec= fread( 'Curve_Fitting_bec.csv')
  sdatabpa= fread( 'Curve_Fitting_bpa.csv')
  sdatatcs= fread( 'Curve_Fitting_tcs.csv')
  sdatami= fread( 'Curve_Fitting_mi.csv')
  sdatamii= fread( 'Curve_Fitting_mii.csv')
  
  sdatabec = sdatabec[xs>med,]
  pvaluebec = mean(sdatabec$PValue)
  sdatabec$treatment_1 = 'BEC'
  sdatabpa = sdatabpa[xs>med,]
  pvaluebpa = mean(sdatabpa$PValue)
  sdatabpa$treatment_1 = 'BPA'
  sdatatcs = sdatatcs[xs>med,]
  pvaluetcs = mean(sdatatcs$PValue)
  sdatatcs$treatment_1 = 'TCS'
  sdatami = sdatami[xs>med,]
  pvaluemi = mean(sdatami$PValue)
  sdatami$treatment_1 = 'MI'
  sdatamii = sdatamii[xs>med,]
  pvaluemii = mean(sdatamii$PValue)
  sdatamii$treatment_1 = 'MII'
  if(pvaluebec ==0){pvaluebec = ' < 0.002'}else{
    pvaluebec = paste('P = ',signif(pvaluebec,2),sep='')
  }
  if(pvaluebpa ==0){pvaluebpa = ' < 0.002'}else{
    pvaluebpa = paste('P = ',signif(pvaluebpa,2),sep='')
  }
  if(pvaluetcs ==0){pvaluetcs = ' < 0.002'}else{
    pvaluetcs = paste('P = ',signif(pvaluetcs,2),sep='')
  }
  if(pvaluemi ==0){pvaluemi = ' < 0.002'}else{
    pvaluemi = paste('P = ',signif(pvaluemi,2),sep='')
  }
  if(pvaluemii ==0){pvaluemii = ' < 0.002'}else{
    pvaluemii = paste('P = ',signif(pvaluemii,2),sep='')
  }
  sdata = rbind(sdatabec,sdatabpa,sdatatcs,sdatami,sdatamii)
  lab_data = data.frame(treatment_1 = c('BEC','BPA','TCS','MI','MII'),
                        lab = paste('\t\t\t\t m =', c(slopebec,slopebpa,slopetcs,slopemi,slopemii),' ' , c(pvaluebec,pvaluebpa,pvaluetcs,pvaluemi,pvaluemii)) ,
                        Overlap = c(med+0.03,med+0.03,med+0.03,med+0.03,med+0.03),Dissimilarity = 0.9)
  dat$treatment_1 = factor(dat$treatment_1,levels=c('BEC','BPA','TCS','MI','MII'))
  lab_data$treatment_1 = factor(lab_data$treatment_1,levels=c('BEC','BPA','TCS','MI','MII'))
  sdata$treatment_1 = factor(sdata$treatment_1,levels=c('BEC','BPA','TCS','MI','MII'))
  p<- ggplot(dat[Overlap>med]) +
    geom_point(data=dat[Overlap>med],aes(x= Overlap,y = Dissimilarity,fill=treatment_1),size=2,color='black',shape=21) +scale_fill_manual(values =alpha(c('black',brewer.pal(5, "Set1")[c(1,2,4,5)]),0.5))+
    theme_bw() + theme_gray(base_family = "Times")+
    geom_line(data=sdata,aes(x=xs,y=mean_lm),col='grey',size=1) +
    geom_ribbon(data=sdata,aes(x=xs,ymin=LCI_lm,ymax=UCI_lm),alpha=0.5,fill='grey') +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.8)) +# theme(axis.line = element_line(size=1)) +
    geom_text(data = lab_data,size=3,family = 'Times',aes(x=Overlap,y=0.76,label=lab)) +
    facet_wrap(~treatment_1,ncol=5,labeller = as_labeller(lab_tr))  + scale_x_continuous(limits=c(0.91,1),breaks=c(0.91,1)) +
    theme(axis.title.y =element_text(face='bold',size=12), axis.title.x=element_text(face='bold',size=12),axis.text.y=element_text(size=10),axis.text.x=element_text(size=10),strip.text = element_text(face='bold',size = 10),legend.title = element_text(size = 12), plot.title = element_text(face='bold',size = 12,hjust = 0.5,vjust = -2), 
          plot.margin = margin(5.5, 5.5, 5.5, 8.5, "pt"),legend.text = element_text(size = 8),legend.key.size = unit(0.2,"cm"),legend.position = "bottom",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=FALSE)+
    theme(panel.spacing = unit(0.8, "lines"))
    
  return(p)
}


doc_df = fread('Dissimilarity_Overlap.csv')
doc_df = doc_df[Dissimilarity>0,]
doc_df2 = doc_df[Same_Environment ==TRUE,]

pB <- plot_lowess(doc_df2,'Curve_Fitting_SameEnv.csv',median(doc_df2$Overlap))
pBInset <- plot_olm(doc_df2,'Curve_Fitting_SameEnv.csv',median(doc_df2$Overlap))+  
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.text =  element_blank(),axis.title = element_blank())
fig3b <- plot_olm_cs(doc_df2,median(doc_df2$Overlap))

#to plot Fig.3, combine fig3a and fig3b, see 'micropollutant_mc_R.R'
  

