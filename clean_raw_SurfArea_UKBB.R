# calculate SA for the 11 regions based on fragmented Desikan atlas
## based on Yash Patel, modified by Zhijie Liao
library(tidyverse)
library(psych)
setwd("~/Desktop/CNV_brain")

## functions to process data:
## clean data: clean outliers defined as +/-1.5 IQR
isnt_out_rIQR <- function(x, thres = 3, na.rm = TRUE) { # QC checking function ... 
  x[which(x < 0.01)] = NA # remove values below ZERO... not possible, mathematically, likely failures at that ROI (possible with finer parc)
  # all areal measures are positive (curv is integrated and rectified through FS)
  #print(table(is.na(x)))
  getrid <- which(x %in% boxplot.stats(x)$out) # basic IQR  # using this for final... most conservative; have worked with adjBoxstats prev. 
  x[getrid] = NA 
  #print(table(is.na(x)))
  return(x)
}
### merge to 11 regions:
## the parcellations in the 11 regions:
region22.list<-list(
  LH.STC = c("LH_parc_144","LH_parc_117", "LH_parc_124", "LH_parc_108"),
  LH.A1C = c("LH_parc_93","LH_parc_111","LH_parc_131"),
  LH.ITC = c("LH_parc_15","LH_parc_16", "LH_parc_7"),
  LH.IPC = c("LH_parc_222","LH_parc_214", "LH_parc_198", "LH_parc_220"),
  LH.S1C = c("LH_parc_275","LH_parc_267","LH_parc_287", "LH_parc_298"),
  LH.M1C = c("LH_parc_308","LH_parc_297","LH_parc_296", "LH_parc_309"),
  LH.VFC = c("LH_parc_239","LH_parc_215","LH_parc_210", "LH_parc_185"),
  LH.DFC = c("LH_parc_233","LH_parc_203","LH_parc_226", "LH_parc_259"),
  LH.OFC = c("LH_parc_99","LH_parc_105", "LH_parc_87", "LH_parc_74"),
  LH.V1C = c("LH_parc_68","LH_parc_73", "LH_parc_57", "LH_parc_70"),
  LH.MFC = c("LH_parc_142","LH_parc_123", "LH_parc_92"),
  
  RH.STC = c("RH_parc_91","RH_parc_73", "RH_parc_89", "RH_parc_113"),
  RH.A1C = c("RH_parc_79","RH_parc_100"),
  RH.ITC = c("RH_parc_10","RH_parc_6", "RH_parc_2"),
  RH.IPC = c("RH_parc_194","RH_parc_170", "RH_parc_187", "RH_parc_172", "RH_parc_162"),
  RH.S1C = c("RH_parc_280","RH_parc_255","RH_parc_266", "RH_parc_259", "RH_parc_237"),
  RH.M1C = c("RH_parc_287","RH_parc_271","RH_parc_296", "RH_parc_277"),
  RH.VFC = c("RH_parc_223","RH_parc_183", "RH_parc_196", "RH_parc_158"),
  RH.DFC = c("RH_parc_241","RH_parc_262","RH_parc_230", "RH_parc_269"),
  RH.OFC = c("RH_parc_86","RH_parc_96", "RH_parc_116", "RH_parc_123"),
  RH.V1C = c("RH_parc_66","RH_parc_83", "RH_parc_87", "RH_parc_71"),
  RH.MFC = c("RH_parc_178","RH_parc_154", "RH_parc_124") )
rois.22<-names(region22.list) ## 22 regions list

get.22roi<-function(clean.df, rois, roi.parc.list){
  df.22roi<-tibble()
  for (i in rois) {
    r<-roi.parc.list[[i]]
    print(r)
    df<-clean.df%>%filter(ROI %in% r)   %>% 
      group_by(sub_id) %>% summarise_at(c("SA","NVtx","Vol","meancurv", "gauscurv", "fold", "curvind", "totalArea","meanArea","ArealDens") ,mean, na.rm=TRUE) %>%
      add_column(ROI = i)
    #print(head(df))
    df.22roi<-rbind(df.22roi, df)
  }
  df.22roi
}


#########################
## ukbb, raw data ######
ukbb.lh<-rbind(read.table("data/ukbb/lh.500vtx.ctxmeasures_chunk1_run2.txt", header = F),
               read.table("data/ukbb/lh.500vtx.ctxmeasures_chunk1_run3.txt", header = F),
               read.table("data/ukbb/lh.500vtx.ctxmeasures_chunk2_run1.txt", header = F),
               read.table("data/ukbb/lh.500vtx.ctxmeasures_chunk2_run2.txt", header = F),
               read.table("data/ukbb/lh.500vtx.ctxmeasures_chunk3.txt", header = F),
               read.table("data/ukbb/lh.500vtx.ctxmeasures_chunk4.txt", header = F) )%>%
  mutate(V11=paste("LH_", V11, sep = ""))
head(ukbb.lh)
length(unique(ukbb.lh$V1))
ukbb.rh<-rbind(read.table("data/ukbb/rh.500vtx.ctxmeasures_chunk1_run2.txt", header = F),
               read.table("data/ukbb/rh.500vtx.ctxmeasures_chunk1_run3.txt", header = F),
               read.table("data/ukbb/rh.500vtx.ctxmeasures_chunk2_run1.txt", header = F),
               read.table("data/ukbb/rh.500vtx.ctxmeasures_chunk2_run2.txt", header = F),
               read.table("data/ukbb/rh.500vtx.ctxmeasures_chunk3.txt", header = F),
               read.table("data/ukbb/rh.500vtx.ctxmeasures_chunk4.txt", header = F) )%>%
  mutate(V11=paste("RH_", V11, sep = ""))


ukbb.all<-rbind(ukbb.lh, ukbb.rh)
names(ukbb.all) <- c("sub_id", "NVtx", "SA", "Vol", "Th", "Th.sd", "meancurv", "gauscurv", "fold", "curvind", "ROI")
ukbb.all<- ukbb.all %>% separate(col="sub_id",into = c("sub_id","a2","a3","a4"))%>%select(-a2, -a3, -a4)%>%
  group_by(sub_id) %>% mutate(totalArea = sum(SA), meanArea= mean(SA))%>%
  group_by(sub_id, ROI) %>% mutate(ArealDens = SA/NVtx)
head(ukbb.all)
##clean outliers, ukbb
ukbb.clean <- ukbb.all %>% group_by(ROI) %>% mutate_at(c("SA","Vol","meancurv", "gauscurv", "fold", "curvind","totalArea","meanArea", "ArealDens") ,isnt_out_rIQR)
table(is.na(ukbb.clean))# 

## get the data of 11 regions, ukbb 
ukbb.22roi<-get.22roi(ukbb.clean, rois.22, region22.list)
head(ukbb.22roi); table(is.na(ukbb.22roi$SA))

# creating left right hemi average data, ukbb
ukbb.11roi<-ukbb.22roi%>%separate(col = ROI, into = c("hemi","roi"))%>%
  group_by(sub_id, roi) %>% summarise_at(c("SA","NVtx","Vol","meancurv", "gauscurv", "fold", "curvind","totalArea","ArealDens") ,sum, na.rm=F)%>%
  mutate(totalArea=totalArea/2)## total area shouldn't be sum of left&right
head(ukbb.11roi)
table(is.na(ukbb.11roi))# 
summary.11roi.ukbb<-ukbb.11roi%>%group_by(roi)%>%summarise_at(c("SA","NVtx","Vol","meancurv", "gauscurv", "fold", "curvind","totalArea","ArealDens") ,mean, na.rm=TRUE)
ggplot(ukbb.11roi, aes(x=SA, color=roi))+
  geom_density()+theme_classic()+ggtitle("UKBB")


###







