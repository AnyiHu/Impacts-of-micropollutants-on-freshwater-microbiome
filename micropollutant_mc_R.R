#this R script is for the analysis and visualization of manuscript 'Ecological drivers of microbial succession following 1 micropollutant exposure'

##set working directory
setwd("/micropollutant_mc2021")

##set font family
windowsFonts(Times=windowsFont("Times New Roman"))

##load packages
library(ggplot2)
library(ggsignif)
library(reshape2)
library(cowplot)
library(data.table)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)

##import data
#data1
pair16s<-read.csv('pair16s.csv')#Bray Curtis distance between treatment and control
pair16s$Days2<-pair16s$Days
pair16s$Days<-as.factor(pair16s$Days)
pair16s2<-melt(pair16s,id=c('Days','Days2'))

#data2
removal<-read.csv('removal.csv')#removal of micropollutants
removal$chemical<-factor(removal$chemical,levels = c('BPA','TCS','BPS','TCC'))
removal$treatment<-factor(removal$treatment,levels = c('BPA','TCS','MI','MII'))
removal$oder<-factor(removal$oder)

##Fig.2 removal
fig2a<-ggplot(removal[which(removal$treatment%in%c('BPA','TCS')),],aes(x=chemical,y=removal,fill=treatment))+
  geom_violin()+ geom_point(shape=1,size=2)+ 
  scale_fill_manual(values =alpha(brewer.pal(5, "Set1")[c(1,2)],0.5))+
  theme_gray(base_family = "Times")+ylab('Removal rate (%)')+
  theme(axis.title.y =element_text(face='bold',size=12), axis.title.x=element_blank(),axis.text.y=element_text(size=10),axis.text.x=element_text(size=10),strip.text = element_text(face='bold',size = 10),legend.title = element_text(size = 12), plot.title = element_text(face='bold',size = 12,hjust = 0.5,vjust = -2), 
        legend.text = element_text(size = 10),legend.key.size = unit(0.4,"cm"),legend.position = "right",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=FALSE)
fig2b<-ggplot(removal[which(!removal$treatment%in%c('BPA','TCS')),],aes(x=oder,y=removal,fill=treatment))+
  geom_violin()+ geom_point(shape=1,size=2)+ 
  scale_fill_manual(values =alpha(brewer.pal(5, "Set1")[c(4,5)],0.5))+
  theme_gray(base_family = "Times")+ylab(' ')+geom_vline(xintercept=2.5)+
  scale_x_discrete(breaks = c('c', 'd', 'e', 'f', 'g','h'),labels = c('BPA', 'TCS', 'BPA', 'TCS', 'BPS','TCC'))+
  theme(axis.title.y =element_text(face='bold',size=12), axis.title.x=element_blank(),axis.text.y=element_text(size=10),axis.text.x=element_text(size=10),strip.text = element_text(face='bold',size = 10),legend.title = element_text(size = 12), plot.title = element_text(face='bold',size = 12,hjust = 0.5,vjust = -2), 
        legend.text = element_text(size = 10),legend.key.size = unit(0.4,"cm"),legend.position = "right",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=FALSE)

fig2<-plot_grid(fig2a, fig2b,labels=c('A','B'), label_fontfamily = 'Times',align = "h", ncol = 2,rel_widths = c(1:3))

##Fig.3
#fig3a Bray Curtis distance between treatment and control over time
compaired<-list(c('5','10'),c('5','15'),c('5','20'),c('5','25'),c('5','30'),c('5','35'))
#label_x<-c('5\n I1','10\n I2','15\n I3','20\n I4','25\n I5','30\n I6','35\n I7')
label_x<-c('B1','B2','B3','B4','B5','B6','B7')
fig3a<-ggplot(pair16s2,aes(x=Days,y=value))+geom_boxplot(color='black')+geom_point(aes(fill=variable), shape=21, color='black',size=2)+
  facet_wrap(.~variable,nrow=1,labeller = as_labeller(c('BPA'='BPA vs. Con.','TCS'='TCS vs. Con.','MI'='MI vs. Con.','MII'='MII vs. Con.')))+
  theme_gray(base_family = "Times")+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = wilcox.test,color='black',textsize = 2,family = 'Times')+
  scale_x_discrete(labels= label_x)+scale_fill_manual(values =alpha(brewer.pal(5, "Set1")[c(1,2,4,5)],0.5))+
  theme(axis.title.y =element_text(face='bold',size=12), axis.title.x=element_text(face='bold',size=12),axis.text.y=element_text(size=10),axis.text.x=element_text(size=10),strip.text = element_text(face='bold',size = 10),legend.title = element_text(size = 12), plot.title = element_text(face='bold',size = 12,hjust = 0.5,vjust = -2), 
        plot.margin = margin(5.5, 8.5, 5.5, 5.5, "pt"),legend.text = element_text(size = 8),legend.key.size = unit(0.2,"cm"),legend.position = "bottom",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(y='Bray-Curtis dissimilarity',x='Batches')+ylim(c(0.25,1.18))+
  guides(fill=FALSE)
#figure3 plot fig3a with fig3b ('5.fig3b_lis.R')
fig3<-plot_grid(fig3a, fig3b,labels=c('A','B'), label_fontfamily = 'Times',align = "h", ncol = 1)

##Fig4 see 'fig4_lis.R'

##Fig5 ecological proccess null test
#the R script of null test is not provided here, please refer to:
#Stegen, J. C., Lin, X., Fredrickson, J. K., & Konopka, A. E. (2015). Estimating and mapping ecological processes influencing microbial community assembly. Frontiers in microbiology, 6, 370.
null<-read.csv('null test_control vs. treatment.csv')
null<-melt(null,id=c('Treatment','Ecological.process'))
null$variable<-factor(null$variable,levels = c('P1','P2','P3'))
null$Treatment<-factor(null$Treatment,levels = c('BPA','TCS','MI','MII'))
null$Ecological.process<-factor(null$Ecological.process,levels = c('Homo_selection','Vari_selection','Dispersal_homogen','Dispersal_limit','Undominated'))
label_x2<-c('P1'='Phase 1','P2'='Phase 2','P3'='Phase 3')
fig5<-ggplot(data=null[which(null$Ecological.process!='Dispersal_homogen'),],aes(x=variable,y=value*100,fill=Ecological.process,group=Ecological.process))+
  theme_gray(base_family = "Times")+
  geom_bar(position = position_stack(reverse = TRUE),colour=NA,stat="identity")+
  facet_wrap(.~Treatment,nrow=2,labeller = as_labeller(c('BPA'='BPA vs. Con.','TCS'='TCS vs. Con.','MI'='MI vs. Con.','MII'='MII vs. Con.')))+
  xlab("Phase")+ylab('Relative proportion (%)')+scale_x_discrete(labels= label_x2)+
  scale_fill_manual(values=brewer.pal(6,"Dark2")[c(1:3,6)],name="Ecological process",labels=c('Homogenizing selection','Variable selection','Dispersal limitation','Ecological drift'))+
  theme(axis.title.y =element_text(face='bold',size=12), axis.title.x=element_blank(),axis.text.y=element_text(size=10),axis.text.x=element_text(size=10),strip.text = element_text(face='bold',size = 10),legend.title = element_text(size = 12), plot.title = element_text(face='bold',size = 12,hjust = 0.5,vjust = -2), 
        legend.text = element_text(size = 10),legend.key.size = unit(0.4,"cm"),legend.position = "right",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())


##Fig6a heatmap see 'ecological group.R'
