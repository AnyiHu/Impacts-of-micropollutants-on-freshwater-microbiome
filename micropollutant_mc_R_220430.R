#this R script is for the analysis and visualization of manuscript 'Repeated introduction of micropollutants enhances microbial succession despite stable degradation patterns'

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
library(phyloseq)
library(vegan)

##import data
#data1
pair16s<-read.csv('pair16s.csv')#Bray Curtis distance between treatment and control
pair16s$Days2<-pair16s$Days
pair16s$Days<-as.factor(pair16s$Days)
pair16s2<-melt(pair16s,id=c('Days','Days2'))

#data2
removal<-read.csv('PPCP removal.csv')#removal of micropollutants
removal$treatment<-factor(removal$treatment,levels = c('BPA','TCS','MI','MII'))
removal$micropollutant<-factor(removal$micropollutant,levels = c('BPA','TCS','BPS','TCC'))

#data3
ls_qiime1 = readRDS("ls_qiime1.rds")
#etract nmds1&2
li_nmds <- ordinate(physeq = ls_qiime1, method = "NMDS", distance = "bray")
li_nmds2<-data.frame(NMDS1=li_nmds$points[,1], NMDS2=li_nmds$points[,2], ls_group)
li_nmds2$Treatment = factor(li_nmds2$Treatment,levels=c("Origin","BEC","BPA","TCS","MI","MII"))
li_nmds2$Time = factor(li_nmds2$Time)
li_nmds2$Batch = factor(li_nmds2$Batch)


##Fig.2 removal over batches
fig2<-ggplot(removal, aes(x=batch, y=removal, color=micropollutant,group=micropollutant))+facet_grid(rows = vars(treatment))+
  geom_smooth(method = 'loess',fill='grey80')+
  geom_jitter(width = 0.2,size=2)+
  scale_color_manual(values =alpha(brewer.pal(5, "Set1")[c(1,2,4,5)],0.5))+
  xlab('Batches')+ylab('Removal (%)')+ 
  scale_x_discrete(breaks = c('I1', 'I2', 'I3', 'I4', 'I5','I6','I7'),labels = c('B1', 'B2', 'B3', 'B4', 'B5','B6','B7'))+
  theme_gray(base_family = "Times")+
  theme(axis.title.y =element_text(face='bold',size=12), axis.title.x=element_text(face='bold',size=12),axis.text.y=element_text(size=10),axis.text.x=element_text(size=10),strip.text = element_text(face='bold',size = 10),legend.title = element_text(size = 12), plot.title = element_text(face='bold',size = 12,hjust = 0.5,vjust = -2), 
        legend.text = element_text(size = 10),legend.key.size = unit(0.4,"cm"),legend.position = "right",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())


##Fig.3 NMDS plot
#nmds using extracted data:li_nmds2
cols <- inlmisc::GetTolColors(n = 10)
labels_Rsq1 <- paste0("Treatment + Time")
labels_Rsq2 <- paste0("PERMANOVA: P < 0.001")
labels_Rsq3 <- paste0("Stress: 0.147")

fig3 <-ggplot() + 
  geom_point(data=li_nmds2, aes(NMDS1, NMDS2, shape=Treatment, fill=Time,color=Time),size=4) + 
  scale_fill_manual(name = "Batches",
                    values  = cols[1:8],  #legend
                    breaks = c(0, 5, 10, 15, 20, 25, 30, 35), 
                    labels = c("Initial", "B1", "B2", "B3", "B4", "B5", "B6", "B7"))  + 
  scale_shape_manual(values = c(20, 21, 22, 23, 24, 25),                              
                     breaks = c("Origin","BEC","BPA","TCS","MI", "MII"), 
                     labels = c("Initial", "Control","BPA","TCS","MI", "MII")) +
  scale_color_manual(values  = c('grey', rep('white',7)))+ #legend
  theme_bw()+ 
  theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())+ 
  theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) + 
  theme(axis.title.y=element_text(colour = 'black',size = 14,family='Times')) + 
  theme(axis.title.x=element_text(colour = 'black',size = 14, family='Times')) + 
  theme(axis.text.y=element_text(colour = 'black', size = 14,family='Times' )) + 
  theme(axis.text.x=element_text(colour = 'black', size = 14, family='Times')) + 
  theme(legend.text=element_text(size=13,family="Times")) + 
  theme(legend.key=element_rect(fill="transparent", colour=NA)) + 
  theme(legend.background=element_rect(fill=NA)) + 
  theme(legend.title=element_text(color="black", size=14, face='bold', family='Times')) + 
  guides(fill=guide_legend(override.aes=list(shape=22, size=5,color='white'))) + 
  guides(shape=guide_legend(override.aes=list(size=4))) + 
  guides(color=FALSE)+
  theme(legend.key.size=unit("0.5","cm")) 
fig3 = fig3 + annotate("text", x =-1.55, y=0.74, size=5, 
                   label=labels_Rsq1, 
                   color="black",parse=TRUE,family="Times") +
  annotate("text", x = -1.3, y=0.66, size=5, 
           label=labels_Rsq2, 
           color="black",parse=TRUE,family="Times")+
  annotate("text", x = -1.7, y=-0.5, size=5, 
           label=labels_Rsq3, 
           color="black",parse=TRUE,family="Times")

#fit with removal
#organize removal
rem<-dcast(removal,treatment+batch+replicate~micropollutant,value.var = 'removal')
#remove that have four NA or imcomplete data
rem<-rem[-c(34,55,61),]
# NA set to 0
rem[is.na(rem)]<-0
#add rem to otu data
write.csv(rem,'rem.csv')

#extract metadata, check the order of samples
ls_group<-data.frame(Time=ls_qiime1@sam_data[["Time3"]],Batch=ls_qiime1@sam_data[["Time"]], Treatment=ls_qiime1@sam_data[["Treatment"]])

#(order manually rem to sample order of ls_group)
x<-read.csv('removal for nmds.csv')
evn<-x[,14:17]

#envfit
fit = envfit(li_nmds, evn, perm = 999,na.rm = TRUE)
arrowdata<-as.data.frame(scores(fit,display="bp"))
arrowdata$micropollutant<-row.names(arrowdata)
arrowdata1<-cbind(x[1:4,1:13],arrowdata)

fig3=fig3+geom_segment(data=arrowdata1,aes(x=0,y=0,xend=NMDS1,yend=NMDS2),colour="black", size=1,linetype=1, arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(data=arrowdata1,mapping=aes(x = NMDS1+0.2,y=NMDS2,label=micropollutant),parse = TRUE,inherit.aes = FALSE,size = 5, family='Times')


##Fig.4
#fig4a Bray Curtis distance between treatment and control over time
compaired<-list(c('5','10'),c('5','15'),c('5','20'),c('5','25'),c('5','30'),c('5','35'))
#label_x<-c('5\n I1','10\n I2','15\n I3','20\n I4','25\n I5','30\n I6','35\n I7')
label_x<-c('B1','B2','B3','B4','B5','B6','B7')
fig4a<-ggplot(pair16s2,aes(x=Days,y=value))+geom_boxplot(color='black')+geom_point(aes(fill=variable), shape=21, color='black',size=2)+
  facet_wrap(.~variable,nrow=1,labeller = as_labeller(c('BPA'='BPA vs. Con.','TCS'='TCS vs. Con.','MI'='MI vs. Con.','MII'='MII vs. Con.')))+
  theme_gray(base_family = "Times")+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = wilcox.test,color='black',textsize = 2,family = 'Times')+
  scale_x_discrete(labels= label_x)+scale_fill_manual(values =alpha(brewer.pal(5, "Set1")[c(1,2,4,5)],0.5))+
  theme(axis.title.y =element_text(face='bold',size=12), axis.title.x=element_text(face='bold',size=12),axis.text.y=element_text(size=10),axis.text.x=element_text(size=10),strip.text = element_text(face='bold',size = 10),legend.title = element_text(size = 12), plot.title = element_text(face='bold',size = 12,hjust = 0.5,vjust = -2), 
        plot.margin = margin(5.5, 8.5, 5.5, 5.5, "pt"),legend.text = element_text(size = 8),legend.key.size = unit(0.2,"cm"),legend.position = "bottom",legend.direction = "vertical",panel.background = element_rect(color = 'black',fill='white'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(y='Bray-Curtis dissimilarity',x='Batches')+ylim(c(0.25,1.18))+
  guides(fill=FALSE)
#figure4 plot fig4a with fig4b ('5.fig4b_lis.R')
fig4<-plot_grid(fig4a, fig4b,labels=c('A','B'), label_fontfamily = 'Times',align = "h", ncol = 1)


##Fig5 heatmap see 'ecological group.R'


##Fig6 ecological proccess null test
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



