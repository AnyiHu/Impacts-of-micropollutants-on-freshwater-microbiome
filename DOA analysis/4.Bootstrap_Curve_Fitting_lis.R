#this manuscript is modified from https://doi.org/10.5281/zenodo.3817698. 
#Vila, J. C., Liu, Y. Y., & Sanchez, A. (2020). Dissimilarity-Overlap analysis of replicate enrichment communities. The ISME Journal, 14(10), 2505-2513.

setwd("/micropollutant_mc2021/DOC analysis")
library(data.table)

calc_lowess<- function(data){
  #fit LOWESS To data and return predicted data.points
  xs = seq(0,1,length=501)
  m <- loess(data$Dissimilarity ~ data$Overlap,span=0.2,family="symmetric",iterations=5,statistics='none')
  y = predict(m,xs)
  return(y)
}

calc_lm<- function(data,med){
  #fit lm to data (above median overlap) and return predicted datapoints.
  xs = seq(0,1,length=501)
  data = data[Overlap>med]
  m <- lm(Dissimilarity ~ Overlap,data)
  #calculate a set of predicted value on row
  y = xs*m$coefficients[2] + m$coefficients[1]
  return(y)
}

calc_lci <- function(sdat,conf){
  return(apply(sdat,2,function(x)  quantile(na.omit(x),(1-conf)/2)))
}

calc_uci <- function(sdat,conf){
  return(apply(sdat,2,function(x)  quantile(na.omit(x),1- (1-conf)/2)))
}

frac_pos_slope <- function(sdat){
  #if m>0, sdat[,ncol(sdat)] > sdat [,ncol(sdat)-1]
  return(sum(sdat[,ncol(sdat)] > sdat [,ncol(sdat)-1])/nrow(sdat)) 
}

calc_frac <-function(data){
  return(nrow(data[Overlap>0.5 & Dissimilarity<sqrt(log(2))/2])/nrow(data))
}

is_pos_slope =function(data,med){
  y = calc_lm(data,med)
  val = y[length(y)] > y[length(y)-1]
  return(val)
}

save_summary <-function(olm,lowess,fn){
  xs = seq(0,1,length=501)
  sdata= data.frame(xs,mean_lm=colMeans(olm,na.rm=TRUE),
                    mean_lowess=colMeans(lowess,na.rm=TRUE),
                    UCI_lm = calc_uci(olm,0.95),
                    UCI_lowess = calc_uci(lowess,0.95),
                    LCI_lm = calc_lci(olm,0.95),
                    LCI_lowess = calc_lci(lowess,0.95),
                    PValue = frac_pos_slope(olm))
  fwrite(sdata,fn)
  
  return()
}


generate_curves <- function(dat,n,null){
  # FITS Lowess and OLM for every bootstrap realization for all of the different subests plotted in the paper (each chunk
  # IS fitted with respect to a one median value)
  #Same Environment broken down by Csources pair  and by same inouculum for figure 1 and S1
  SameEnv_olm = c()
  SameEnv_lowess = c()
  bec_olm <- c()
  bec_lowess<- c()
  bpa_olm<- c()
  bpa_lowess<- c()
  tcs_olm<- c()
  tcs_lowess<- c()
  mi_olm<- c()
  mi_lowess<- c()
  mii_olm<- c()
  mii_lowess<- c()
  
  
  #Different Environment broken down by Csources pair for figure s2 and s3
  DiffEnv_olm = c()
  DiffEnv_lowess = c()
  bec_bpa_olm <- c()
  bec_bpa_lowess<- c()
  bec_tcs_olm<- c()
  bec_tcs_lowess<- c()
  bec_mi_olm<- c()
  bec_mi_lowess<- c()
  bec_mii_olm<- c()
  bec_mii_lowess<- c()
  bpa_tcs_olm<- c()
  bpa_tcs_lowess<- c()
  bpa_mi_olm<- c()
  bpa_mi_lowess<- c()
  bpa_mii_olm<- c()
  bpa_mii_lowess<- c()
  tcs_mi_olm<- c()
  tcs_mi_lowess<- c()
  tcs_mii_olm<- c()
  tcs_mii_lowess<- c()
  mi_mii_olm<- c()
  mi_mii_lowess<- c()
  
  
  #Break down by different combinations of same environment and same inoculum for figure s10
  All_olm<- c() #All
  All_lowess<- c() #All
 

  for(j in 1:n){
    print(j)
    #Data for all pairs
    All = dat[Run==j,]
    
    #First data for same environment
    SameEnvdf = All[Same_Environment == TRUE]
    becdf <- SameEnvdf[treatment_1 == 'BEC',]
    bpadf <- SameEnvdf[treatment_1 == 'BPA',]
    tcsdf <- SameEnvdf[treatment_1 == 'TCS',]
    midf <- SameEnvdf[treatment_1 == 'MI',]
    miidf <- SameEnvdf[treatment_1 == 'MII',]
    
    #next data for different environment
    DiffEnvdf = All[Same_Environment == FALSE]
    becbpadf <- DiffEnvdf[treatment_1 == 'BPA' & treatment_2 == 'BEC',]
    bectcsdf <- DiffEnvdf[treatment_1 == 'TCS' & treatment_2 == 'BEC',]
    becmidf <- DiffEnvdf[treatment_1 == 'MI' & treatment_2 == 'BEC',]
    becmiidf <- DiffEnvdf[treatment_1 == 'MII' & treatment_2 == 'BEC',]
    bpatcsdf <- DiffEnvdf[treatment_1 == 'TCS' & treatment_2 == 'BPA',]
    bpamidf <- DiffEnvdf[treatment_1 == 'MI' & treatment_2 == 'BPA',]
    bpamiidf <- DiffEnvdf[treatment_1 == 'MII' & treatment_2 == 'BPA',]
    tcsmidf <- DiffEnvdf[treatment_1 == 'TCS' & treatment_2 == 'MI',]
    tcsmiidf <- DiffEnvdf[treatment_1 == 'TCS' & treatment_2 == 'MII',]
    mimiidf <- DiffEnvdf[treatment_1 == 'MII' & treatment_2 == 'MI',]
   
    #first fit data for same environment (using median for same environment pairs)
    SameEnv_olm<- rbind(SameEnv_olm,calc_lm(SameEnvdf,median(SameEnvdf$Overlap)))
    SameEnv_lowess<- rbind(SameEnv_lowess,calc_lowess(SameEnvdf))
    bec_olm<- rbind(bec_olm,calc_lm(becdf,median(SameEnvdf$Overlap)))
    bec_lowess<- rbind(bec_lowess,calc_lowess(becdf))
    bpa_olm<- rbind(bpa_olm,calc_lm(bpadf,median(SameEnvdf$Overlap)))
    bpa_lowess<- rbind(bpa_lowess,calc_lowess(bpadf))
    tcs_olm<- rbind(tcs_olm,calc_lm(tcsdf,median(SameEnvdf$Overlap)))
    tcs_lowess<- rbind(tcs_lowess,calc_lowess(tcsdf))
    mi_olm<- rbind(mi_olm,calc_lm(midf,median(SameEnvdf$Overlap)))
    mi_lowess<- rbind(mi_lowess,calc_lowess(midf))
    mii_olm<- rbind(mii_olm,calc_lm(miidf,median(SameEnvdf$Overlap)))
    mii_lowess<- rbind(mii_lowess,calc_lowess(miidf))
    
    #next fit data for different environment (using median for different environment pairs)
    DiffEnv_olm<- rbind(DiffEnv_olm,calc_lm(DiffEnvdf,median(DiffEnvdf$Overlap)))
    DiffEnv_lowess<- rbind(DiffEnv_lowess,calc_lowess(DiffEnvdf))
    bec_bpa_olm<- rbind(bec_bpa_olm,calc_lm(becbpadf,median(DiffEnvdf$Overlap)))
    bec_bpa_lowess<- rbind(bec_bpa_lowess,calc_lowess(becbpadf))
    bec_tcs_olm<- rbind(bec_tcs_olm,calc_lm(bectcsdf,median(DiffEnvdf$Overlap)))
    bec_tcs_lowess<- rbind(bec_tcs_lowess,calc_lowess(bectcsdf))
    bec_mi_olm<- rbind(bec_mi_olm,calc_lm(becmidf,median(DiffEnvdf$Overlap)))
    bec_mi_lowess<- rbind(bec_mi_lowess,calc_lowess(becmidf))
    bec_mii_olm<- rbind(bec_mii_olm,calc_lm(becmiidf,median(DiffEnvdf$Overlap)))
    bec_mii_lowess<- rbind(bec_mii_lowess,calc_lowess(becmiidf))
    bpa_tcs_olm<- rbind(bpa_tcs_olm,calc_lm(bpatcsdf,median(DiffEnvdf$Overlap)))
    bpa_tcs_lowess<- rbind(bpa_tcs_lowess,calc_lowess(bpatcsdf))
    bpa_mi_olm<- rbind(bpa_mi_olm,calc_lm(bpamidf,median(DiffEnvdf$Overlap)))
    bpa_mi_lowess<- rbind(bpa_mi_lowess,calc_lowess(bpamidf))
    bpa_mii_olm<- rbind(bpa_mii_olm,calc_lm(bpamiidf,median(DiffEnvdf$Overlap)))
    bpa_mii_lowess<- rbind(bpa_mii_lowess,calc_lowess(bpamiidf))
    tcs_mi_olm<- rbind(tcs_mi_olm,calc_lm(tcsmidf,median(DiffEnvdf$Overlap)))
    tcs_mi_lowess<- rbind(tcs_mi_lowess,calc_lowess(tcsmidf))
    tcs_mii_olm<- rbind(tcs_mii_olm,calc_lm(tcsmiidf,median(DiffEnvdf$Overlap)))
    tcs_mii_lowess<- rbind(tcs_mii_lowess,calc_lowess(tcsmiidf))
    mi_mii_olm<- rbind(mi_mii_olm,calc_lm(mimiidf,median(DiffEnvdf$Overlap)))
    mi_mii_lowess<- rbind(mi_mii_lowess,calc_lowess(mimiidf))
    
    #Finally different combinations of same environment and same inoculum using all is reference point for figure s10
    All_olm<- rbind(All_olm,calc_lm(All,median(All$Overlap))) #All
    All_lowess<-  rbind(All_lowess,calc_lowess(All)) #All
    
  }
  if(null==0){
    #ALL
    save_summary(All_olm,All_lowess,'Curve_Fitting_All.csv')
    

    #SAME ENV
    save_summary(SameEnv_olm,SameEnv_lowess,'Curve_Fitting_SameEnv.csv')
    save_summary(bec_olm,bec_lowess,'Curve_Fitting_bec.csv')
    save_summary(bpa_olm,bpa_lowess,'Curve_Fitting_bpa.csv')
    save_summary(tcs_olm,tcs_lowess,'Curve_Fitting_tcs.csv')
    save_summary(mi_olm,mi_lowess,'Curve_Fitting_mi.csv')
    save_summary(mii_olm,mii_lowess,'Curve_Fitting_mii.csv')
    
    #DIfferent Env
    save_summary(DiffEnv_olm,DiffEnv_lowess,'Curve_Fitting_DiffEnv.csv')
    save_summary(bec_bpa_olm,bec_bpa_lowess,'Curve_Fitting_bec_bpa.csv')
    save_summary(bec_tcs_olm,bec_tcs_lowess,'Curve_Fitting_bec_tcs.csv')
    save_summary(bec_mi_olm,bec_mi_lowess,'Curve_Fitting_bec_mi.csv')
    save_summary(bec_mii_olm,bec_mii_lowess,'Curve_Fitting_bec_mii.csv')
    save_summary(bpa_tcs_olm,bpa_tcs_lowess,'Curve_Fitting_bpa_tcs.csv')
    save_summary(bpa_mi_olm,bpa_mi_lowess,'Curve_Fitting_bpa_mi.csv')
    save_summary(bpa_mii_olm,bpa_mii_lowess,'Curve_Fitting_bpa_mii.csv')
    save_summary(tcs_mi_olm,tcs_mi_lowess,'Curve_Fitting_tcs_mi.csv')
    save_summary(tcs_mii_olm,tcs_mii_lowess,'Curve_Fitting_tcs_mii.csv')
    save_summary(mi_mii_olm,mi_mii_lowess,'Curve_Fitting_mi_mii.csv')
    
  } else{
    #ALL
    save_summary(All_olm,All_lowess,'Curve_Fitting_All_null.csv')
    
    
    #SAME ENV
    save_summary(SameEnv_olm,SameEnv_lowess,'Curve_Fitting_SameEnv_null.csv')
    save_summary(bec_olm,bec_lowess,'Curve_Fitting_bec_null.csv')
    save_summary(bpa_olm,bpa_lowess,'Curve_Fitting_bpa_null.csv')
    save_summary(tcs_olm,tcs_lowess,'Curve_Fitting_tcs_null.csv')
    save_summary(mi_olm,mi_lowess,'Curve_Fitting_mi_null.csv')
    save_summary(mii_olm,mii_lowess,'Curve_Fitting_mii_null.csv')
    
    #DIfferent Env
    save_summary(DiffEnv_olm,DiffEnv_lowess,'Curve_Fitting_DiffEnv_null.csv')
    save_summary(bec_bpa_olm,bec_bpa_lowess,'Curve_Fitting_bec_bpa_null.csv')
    save_summary(bec_tcs_olm,bec_tcs_lowess,'Curve_Fitting_bec_tcs_null.csv')
    save_summary(bec_mi_olm,bec_mi_lowess,'Curve_Fitting_bec_mi_null.csv')
    save_summary(bec_mii_olm,bec_mii_lowess,'Curve_Fitting_bec_mii_null.csv')
    save_summary(bpa_tcs_olm,bpa_tcs_lowess,'Curve_Fitting_bpa_tcs_null.csv')
    save_summary(bpa_mi_olm,bpa_mi_lowess,'Curve_Fitting_bpa_mi_null.csv')
    save_summary(bpa_mii_olm,bpa_mii_lowess,'Curve_Fitting_bpa_mii_null.csv')
    save_summary(tcs_mi_olm,tcs_mi_lowess,'Curve_Fitting_tcs_mi_null.csv')
    save_summary(tcs_mii_olm,tcs_mii_lowess,'Curve_Fitting_tcs_mii_null.csv')
    save_summary(mi_mii_olm,mi_mii_lowess,'Curve_Fitting_mi_mii_null.csv')

  }
  return()
}


set.seed(1)
doc_df_b=fread('Dissimilarity_Overlap_Bootstrapped.csv')
doc_df_b = doc_df_b[Dissimilarity>0,]
generate_curves(doc_df_b,500,0)
set.seed(1)
doc_df_b=fread('Dissimilarity_Overlap_Bootstrapped_Null.csv')
doc_df_b = doc_df_b[Dissimilarity>0,]
generate_curves(doc_df_b,500,1)
