#this manuscript is modified from https://doi.org/10.5281/zenodo.3817698. 
#Vila, J. C., Liu, Y. Y., & Sanchez, A. (2020). Dissimilarity-Overlap analysis of replicate enrichment communities. The ISME Journal, 14(10), 2505-2513.

setwd("/micropollutant_mc2021/DOC analysis")
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)

windowsFonts(Times=windowsFont("Times New Roman"))
#rm(list=ls())

plot_lowess2<- function(dat,fn,med){
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


plot_olm2<- function(dat,fn,med){
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


lab_tr2<-c('BEC-BPA'='Con. vs. BPA','BEC-TCS'='Con. vs. TCS','BEC-MI'='Con. vs. MI','BEC-MII'='Con. vs. MII','BPA-TCS'='BPA vs. TCS','BPA-MI'='BPA vs. MI','BPA-MII'='BPA vs. MII','TCS-MI'='TCS vs. MI','TCS-MII'='TCS vs. MII','MI-MII'='MI vs. MII')
plot_olm_cs2<- function(dat,med){
  #Plot olm over subset of data with above med overlap gationg by carbon source for figure C and supplementary
  # Fn is file containing bootstrapped data
  #Subset
  dat = dat[Overlap>med]
  becbpa <- dat[treatment_1 == 'BPA' & treatment_2 == 'BEC',]
  becbpa$treatment_1 = 'BEC-BPA'
  bectcs <- dat[treatment_1 == 'TCS' & treatment_2 == 'BEC',]
  bectcs$treatment_1 = 'BEC-TCS'
  becmi <- dat[treatment_1 == 'MI' & treatment_2 == 'BEC',]
  becmi$treatment_1 = 'BEC-MI'
  becmii <- dat[treatment_1 == 'MII' & treatment_2 == 'BEC',]
  becmii$treatment_1 = 'BEC-MII'
  bpatcs <- dat[treatment_1 == 'TCS' & treatment_2 == 'BPA',]
  bpatcs$treatment_1 = 'BPA-TCS'
  bpami <- dat[treatment_1 == 'MI' & treatment_2 == 'BPA',]
  bpami$treatment_1 = 'BPA-MI'
  bpamii <- dat[treatment_1 == 'MII' & treatment_2 == 'BPA',]
  bpamii$treatment_1 = 'BPA-MII'
  tcsmi <- dat[treatment_1 == 'TCS' & treatment_2 == 'MI',]
  tcsmi$treatment_1 = 'TCS-MI'
  tcsmii <- dat[treatment_1 == 'TCS' & treatment_2 == 'MII',]
  tcsmii$treatment_1 = 'TCS-MII'
  mimii <- dat[treatment_1 == 'MII' & treatment_2 == 'MI',]
  mimii$treatment_1 = 'MI-MII'
  #media
  mbecbpa = lm(Dissimilarity ~ Overlap, becbpa)
  mbectcs = lm(Dissimilarity ~ Overlap, bectcs)
  mbecmi = lm(Dissimilarity ~ Overlap, becmi)
  mbecmii = lm(Dissimilarity ~ Overlap, becmii)
  mbpatcs = lm(Dissimilarity ~ Overlap, bpatcs)
  mbpami = lm(Dissimilarity ~ Overlap, bpami)
  mbpamii = lm(Dissimilarity ~ Overlap, bpamii)
  mtcsmi = lm(Dissimilarity ~ Overlap, tcsmi)
  mtcsmii = lm(Dissimilarity ~ Overlap, tcsmii)
  mmimii = lm(Dissimilarity ~ Overlap, mimii)
  
  slopebecbpa =floor(signif(mbecbpa$coefficients[2],2)*100)/100
  slopebectcs =floor(signif(mbectcs$coefficients[2],2)*100)/100
  slopebecmi =floor(signif(mbecmi$coefficients[2],2)*100)/100
  slopebecmii =floor(signif(mbecmii$coefficients[2],2)*100)/100
  slopebpatcs =floor(signif(mbpatcs$coefficients[2],2)*100)/100
  slopebpami =floor(signif(mbpami$coefficients[2],2)*100)/100
  slopebpamii =floor(signif(mbpamii$coefficients[2],2)*100)/100
  slopetcsmi =floor(signif(mtcsmi$coefficients[2],2)*100)/100
  slopetcsmii =floor(signif(mtcsmii$coefficients[2],2)*100)/100
  slopemimii =floor(signif(mmimii$coefficients[2],2)*100)/100
  
  sdatabecbpa= fread( 'Curve_Fitting_bec_bpa.csv')
  sdatabectcs= fread( 'Curve_Fitting_bec_tcs.csv')
  sdatabecmi= fread( 'Curve_Fitting_bec_mi.csv')
  sdatabecmii= fread( 'Curve_Fitting_bec_mii.csv')
  sdatabpatcs= fread( 'Curve_Fitting_bpa_tcs.csv')
  sdatabpami= fread( 'Curve_Fitting_bpa_mi.csv')
  sdatabpamii= fread( 'Curve_Fitting_bpa_mii.csv')
  sdatatcsmi= fread( 'Curve_Fitting_tcs_mi.csv')
  sdatatcsmii= fread( 'Curve_Fitting_tcs_mii.csv')
  sdatamimii= fread( 'Curve_Fitting_mi_mii.csv')
  
  sdatabecbpa = sdatabecbpa[xs>med,]
  pvaluebecbpa = mean(sdatabecbpa$PValue)
  sdatabecbpa$treatment_1 = 'BEC-BPA'
  sdatabectcs = sdatabectcs[xs>med,]
  pvaluebectcs = mean(sdatabectcs$PValue)
  sdatabectcs$treatment_1 = 'BEC-TCS'
  sdatabecmi = sdatabecmi[xs>med,]
  pvaluebecmi = mean(sdatabecmi$PValue)
  sdatabecmi$treatment_1 = 'BEC-MI'
  sdatabecmii = sdatabecmii[xs>med,]
  pvaluebecmii = mean(sdatabecmii$PValue)
  sdatabecmii$treatment_1 = 'BEC-MII'
  sdatabpatcs = sdatabpatcs[xs>med,]
  pvaluebpatcs = mean(sdatabpatcs$PValue)
  sdatabpatcs$treatment_1 = 'BPA-TCS'
  sdatabpami = sdatabpami[xs>med,]
  pvaluebpami = mean(sdatabpami$PValue)
  sdatabpami$treatment_1 = 'BPA-MI'
  sdatabpamii = sdatabpamii[xs>med,]
  pvaluebpamii = mean(sdatabpamii$PValue)
  sdatabpamii$treatment_1 = 'BPA-MII'
  sdatatcsmi = sdatatcsmi[xs>med,]
  pvaluetcsmi = mean(sdatatcsmi$PValue)
  sdatatcsmi$treatment_1 = 'TCS-MI'
  sdatatcsmii = sdatatcsmii[xs>med,]
  pvaluetcsmii = mean(sdatatcsmii$PValue)
  sdatatcsmii$treatment_1 = 'TCS-MII'
  sdatamimii = sdatamimii[xs>med,]
  pvaluemimii = mean(sdatamimii$PValue)
  sdatamimii$treatment_1 = 'MI-MII'
  
  if(pvaluebecbpa ==0){pvaluebecbpa = 'p < 0.002'}else{
    pvaluebecbpa = paste('P = ',signif(pvaluebecbpa,2),sep='')
  }
  if(pvaluebectcs ==0){pvaluebectcs = 'p < 0.002'}else{
    pvaluebectcs = paste('P = ',signif(pvaluebectcs,2),sep='')
  }
  if(pvaluebecmi ==0){pvaluebecmi = 'p < 0.002'}else{
    pvaluebecmi = paste('P = ',signif(pvaluebecmi,2),sep='')
  }
  if(pvaluebecmii ==0){pvaluebecmii = 'p < 0.002'}else{
    pvaluebecmii = paste('P = ',signif(pvaluebecmii,2),sep='')
  }
  if(pvaluebpatcs ==0){pvaluebpatcs = 'p < 0.002'}else{
    pvaluebpatcs = paste('P = ',signif(pvaluebpatcs,2),sep='')
  }
  if(pvaluebpami ==0){pvaluebpami = 'p < 0.002'}else{
    pvaluebpami = paste('P = ',signif(pvaluebpami,2),sep='')
  }
  if(pvaluebpamii ==0){pvaluebpamii = 'p < 0.002'}else{
    pvaluebpamii = paste('P = ',signif(pvaluebpamii,2),sep='')
  }
  if(pvaluetcsmi ==0){pvaluetcsmi = 'p < 0.002'}else{
    pvaluetcsmi = paste('P = ',signif(pvaluetcsmi,2),sep='')
  }
  if(pvaluetcsmii ==0){pvaluetcsmii = 'p < 0.002'}else{
    pvaluetcsmii = paste('P = ',signif(pvaluetcsmii,2),sep='')
  }
  if(pvaluemimii ==0){pvaluemimii = 'p < 0.002'}else{
    pvaluemimii = paste('P = ',signif(pvaluemimii,2),sep='')
  }

  sdata = rbind(sdatabecbpa,sdatabectcs,sdatabecmi,sdatabecmii,sdatabpatcs,sdatabpami,sdatabpamii,sdatatcsmi,sdatatcsmii,sdatamimii)
  lab_data = data.frame(treatment_1 = c('BEC-BPA','BEC-TCS','BEC-MI','BEC-MII','BPA-TCS','BPA-MI','BPA-MII','TCS-MI','TCS-MII','MI-MII'),
                        lab = paste('\t\t\t\t m =', c(slopebecbpa,slopebectcs,slopebecmi,slopebecmii,slopebpatcs,slopebpami,slopebpamii,slopetcsmi,slopetcsmii,slopemimii),' ' , c(pvaluebecbpa,pvaluebectcs,pvaluebecmi,pvaluebecmii,pvaluebpatcs,pvaluebpami,pvaluebpamii,pvaluetcsmi,pvaluetcsmii,pvaluemimii)) ,
                        Overlap = c(med+0.05,med+0.05,med+0.05,med+0.05,med+0.05),Dissimilarity = 0.9)
  plot_dat = rbind(becbpa,bectcs,becmi,becmii,bpatcs,bpami,bpamii,tcsmi,tcsmii,mimii)
  sdata$treatment_1 = factor(sdata$treatment_1,levels=c('BEC-BPA','BEC-TCS','BEC-MI','BEC-MII','BPA-TCS','BPA-MI','BPA-MII','TCS-MI','TCS-MII','MI-MII'))
  lab_data$treatment_1 = factor(lab_data$treatment_1,levels=c('BEC-BPA','BEC-TCS','BEC-MI','BEC-MII','BPA-TCS','BPA-MI','BPA-MII','TCS-MI','TCS-MII','MI-MII'))
  plot_dat$treatment_1 = factor(plot_dat$treatment_1,levels=c('BEC-BPA','BEC-TCS','BEC-MI','BEC-MII','BPA-TCS','BPA-MI','BPA-MII','TCS-MI','TCS-MII','MI-MII'))
  
  p<- ggplot(plot_dat) +
    geom_point(data=plot_dat,aes(x= Overlap,y = Dissimilarity),size=2,fill=alpha('black',0.5),color='black',shape=21) +
    theme_bw() + theme_gray(base_family = "Times")+
    geom_line(data=sdata,aes(x=xs,y=mean_lm),col='grey',size=1) +
    geom_ribbon(data=sdata,aes(x=xs,ymin=LCI_lm,ymax=UCI_lm),alpha=0.5,fill='grey') +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.8)) +# theme(axis.line = element_line(size=1)) +
    geom_text(data = lab_data,size=3,family = 'Times',aes(x=Overlap,y=0.76,label=lab)) +
    facet_wrap(~treatment_1,ncol=5,labeller = as_labeller(lab_tr2))  + scale_x_continuous(limits=c(0.85,1),breaks=c(0.85,1)) +
    theme(axis.title.y =element_text(face='bold',size=12), axis.title.x=element_text(face='bold',size=12),axis.text.y=element_text(size=8),axis.text.x=element_text(size=8),strip.text = element_text(face='bold',size = 10),legend.title = element_text(size = 12), plot.title = element_text(face='bold',size = 12,hjust = 0.5,vjust = -2), 
          legend.text = element_text(size = 8),legend.key.size = unit(0.2,"cm"),legend.position = "bottom",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.spacing = unit(0.5, "lines"))+
    guides(color=FALSE)+
    theme(panel.spacing = unit(0.7, "lines"))
  
  return(p)
}



doc_df = fread('Dissimilarity_Overlap.csv')
doc_df = doc_df[Dissimilarity>0,]
doc_df2 = doc_df[Same_Environment ==FALSE,]

pB2 <- plot_lowess2(doc_df2,'Curve_Fitting_DiffEnv.csv',median(doc_df2$Overlap))
pBInset2 <- plot_olm2(doc_df2,'Curve_Fitting_DiffEnv.csv',median(doc_df2$Overlap))+  
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.text =  element_blank(),axis.title = element_blank())
fig4 <- plot_olm_cs2(doc_df2,median(doc_df2$Overlap))



