#this manuscript is modified from https://doi.org/10.5281/zenodo.3817698. 
#Vila, J. C., Liu, Y. Y., & Sanchez, A. (2020). Dissimilarity-Overlap analysis of replicate enrichment communities. The ISME Journal, 14(10), 2505-2513.

setwd("/micropollutant_mc2021/DOC analysis")
library(data.table)
library(operators)

t_test_csource <- function(dat,j,o){
  dat = dat[Run ==j & Overlap>o]
  becbec <- dat[treatment_1 == 'BEC' & treatment_2 == 'BEC',]
  bpabpa <- dat[treatment_1 == 'BPA' & treatment_2 == 'BPA',]
  tcstcs <- dat[treatment_1 == 'TCS' & treatment_2 == 'TCS',]
  mimi <- dat[treatment_1 == 'MI' & treatment_2 == 'MI',]
  miimii <- dat[treatment_1 == 'MII' & treatment_2 == 'MII',]
  becbpa <- dat[treatment_1 == 'BPA' & treatment_2 == 'BEC',]
  bectcs <- dat[treatment_1 == 'TCS' & treatment_2 == 'BEC',]
  becmi <- dat[treatment_1 == 'MI' & treatment_2 == 'BEC',]
  becmii <- dat[treatment_1 == 'MII' & treatment_2 == 'BEC',]
  bpatcs <- dat[treatment_1 == 'TCS' & treatment_2 == 'BPA',]
  bpami <- dat[treatment_1 == 'MI' & treatment_2 == 'BPA',]
  bpamii <- dat[treatment_1 == 'MII' & treatment_2 == 'BPA',]
  tcsmi <- dat[treatment_1 == 'TCS' & treatment_2 == 'MI',]
  tcsmii <- dat[treatment_1 == 'TCS' & treatment_2 == 'MII',]
  mimii <- dat[treatment_1 == 'MII' & treatment_2 == 'MI',]
  t1  =t.test(bpabpa$Dissimilarity,becbec$Dissimilarity)
  t2 = t.test(bpabpa$Dissimilarity,becbpa$Dissimilarity)
  t3 = t.test(becbec$Dissimilarity,becbpa$Dissimilarity)
  t4  =t.test(tcstcs$Dissimilarity,becbec$Dissimilarity)
  t5 = t.test(tcstcs$Dissimilarity,bectcs$Dissimilarity)
  t6 = t.test(becbec$Dissimilarity,bectcs$Dissimilarity)
  t7  =t.test(mimi$Dissimilarity,becbec$Dissimilarity)
  t8 = t.test(mimi$Dissimilarity,becmi$Dissimilarity)
  t9 = t.test(becbec$Dissimilarity,becmi$Dissimilarity)
  t10  =t.test(miimii$Dissimilarity,becbec$Dissimilarity)
  t11 = t.test(miimii$Dissimilarity,becmii$Dissimilarity)
  t12 = t.test(becbec$Dissimilarity,becmii$Dissimilarity)
  t13  =t.test(bpabpa$Dissimilarity,tcstcs$Dissimilarity)
  t14 = t.test(bpabpa$Dissimilarity,bpatcs$Dissimilarity)
  t15 = t.test(tcstcs$Dissimilarity,bpatcs$Dissimilarity)
  t16  =t.test(bpabpa$Dissimilarity,mimi$Dissimilarity)
  t17 = t.test(bpabpa$Dissimilarity,bpami$Dissimilarity)
  t18 = t.test(mimi$Dissimilarity,bpami$Dissimilarity)
  t19  =t.test(bpabpa$Dissimilarity,miimii$Dissimilarity)
  t20 = t.test(bpabpa$Dissimilarity,bpamii$Dissimilarity)
  t21 = t.test(miimii$Dissimilarity,bpamii$Dissimilarity)
  t22  =t.test(tcstcs$Dissimilarity,mimi$Dissimilarity)
  t23 = t.test(tcstcs$Dissimilarity,tcsmi$Dissimilarity)
  t24 = t.test(mimi$Dissimilarity,tcsmi$Dissimilarity)
  t25  =t.test(tcstcs$Dissimilarity,miimii$Dissimilarity)
  t26 = t.test(tcstcs$Dissimilarity,tcsmii$Dissimilarity)
  t27 = t.test(miimii$Dissimilarity,tcsmii$Dissimilarity)
  t28  =t.test(mimi$Dissimilarity,miimii$Dissimilarity)
  t29 = t.test(mimi$Dissimilarity,mimii$Dissimilarity)
  t30 = t.test(miimii$Dissimilarity,mimii$Dissimilarity)
  
  return(data.frame(Comparison = c('bpabpa-becbec','bpabpa-bpabec','bpabec-becbec',
                                   'tcstcs-becbec','tcstcs-tcsbec','tcsbec-becbec',
                                   'mimi-becbec','mimi-mibec','mibec-becbec',
                                   'miimii-becbec','miimii-miibec','miibec-becbec',
                                   'tcstcs-bpabpa','tcstcs-tcsbpa','tcsbpa-bpabpa',
                                   'mimi-bpabpa','mimi-mibpa','mibpa-bpabpa',
                                   'miimii-bpabpa','miimii-miibpa','miibpa-bpabpa',
                                   'mimi-tcstcs','mimi-mitcs','mitcs-tcstcs',
                                   'miimii-tcstcs','miimii-miitcs','miitcs-tcstcs',
                                   'miimii-mimi','miimii-miimi','miimi-mimi'),
    Threshold =o,
    Run =j,
    t = c(t1$statistic,t2$statistic,t3$statistic,
          t4$statistic,t5$statistic,t6$statistic,
          t7$statistic,t8$statistic,t9$statistic,
          t10$statistic,t11$statistic,t12$statistic,
          t13$statistic,t14$statistic,t15$statistic,
          t16$statistic,t17$statistic,t18$statistic,
          t19$statistic,t20$statistic,t21$statistic,
          t22$statistic,t23$statistic,t24$statistic,
          t25$statistic,t26$statistic,t27$statistic,
          t28$statistic,t29$statistic,t30$statistic),
    p = c(t1$p.value,t2$p.value,t3$p.value,
          t4$p.value,t5$p.value,t6$p.value,
          t7$p.value,t8$p.value,t9$p.value,
          t10$p.value,t11$p.value,t12$p.value,
          t13$p.value,t14$p.value,t15$p.value,
          t16$p.value,t17$p.value,t18$p.value,
          t19$p.value,t20$p.value,t21$p.value,
          t22$p.value,t23$p.value,t24$p.value,
          t25$p.value,t26$p.value,t27$p.value,
          t28$p.value,t29$p.value,t30$p.value)))
}


set.seed(1)
doc_df=fread('Dissimilarity_Overlap.csv')
doc_df_b=fread('Dissimilarity_Overlap_Bootstrapped.csv')
doc_df_b_null=fread('Dissimilarity_Overlap_Bootstrapped_Null.csv')
doc_df = doc_df[Dissimilarity>0,]
doc_df_b = doc_df_b[Dissimilarity>0,]
doc_df_b_null = doc_df_b_null[Dissimilarity>0,]

all_df = data.frame()
all_df_null = data.frame()



for(j in 1:max(doc_df_b$Run)){
  print(j)
  all_df_null = rbind(all_df_null,t_test_csource(doc_df_b_null,j,0.9))
  for(o in c(0.9,0.91,0.92,0.93,0.94,0.95)){
    
    all_df = rbind(all_df,tryCatch(t_test_csource(doc_df_b,j,o),
                                   error=function(e){return(data.frame(Comparison = c('bpabpa-becbec','bpabpa-bpabec','bpabec-becbec',
                                                                                      'tcstcs-becbec','tcstcs-tcsbec','tcsbec-becbec',
                                                                                      'mimi-becbec','mimi-mibec','mibec-becbec',
                                                                                      'miimii-becbec','miimii-miibec','miibec-becbec',
                                                                                      'tcstcs-bpabpa','tcstcs-tcsbpa','tcsbpa-bpabpa',
                                                                                      'mimi-bpabpa','mimi-mibpa','mibpa-bpabpa',
                                                                                      'miimii-bpabpa','miimii-miibpa','miibpa-bpabpa',
                                                                                      'mimi-tcstcs','mimi-mitcs','mitcs-tcstcs',
                                                                                      'miimii-tcstcs','miimii-miitcs','miitcs-tcstcs',
                                                                                      'miimii-mimi','miimii-miimi','miimi-mimi'),
                                                                                                    Threshold =o,
                                                                                                    Run =j,
                                                                                                    t = NA,
                                                                                                    p = NA))}))
    
  }
}

fwrite(all_df,'TTest_CSource.csv')
fwrite(all_df_null,'TTest_CSource_Null.csv')



