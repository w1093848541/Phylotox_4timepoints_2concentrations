#BiocManager::install("ggsignif")
#BiocManager::install("patchwork")
#a
suppressPackageStartupMessages({
  # Bioconductor packages
  library(structToolbox)
  library(pmp)
  library(ropls)
  library(BiocFileCache)
  
  # CRAN libraries
  library(ggplot2)
  library(gridExtra)
  library(cowplot)
  library(openxlsx)
  library(ggthemes)
  library(RColorBrewer)
  library(ropls)
  library(dplyr)
  library(mixOmics)
})
setwd("C:/Users/10938/Desktop/Jacobson Lab/2021.10-/4 time points")
MATLABresult = read.xlsx("Combined_withQC_4tp_all.xlsx",sheet = 1,colNames= FALSE,rowNames = FALSE);
MassList = MATLABresult$X1;
MassList = MassList[-1];
Mass_list = as.character(MassList);
PhylotoxData = MTBLS79
NumSamples = as.numeric(dim(MATLABresult)[2])-1
NumMetabolites = as.numeric(dim(MATLABresult)[1])-1
nSamples = 1
RawMatrix = as.matrix(MATLABresult)
Type_names = vector();
while (nSamples <= NumSamples) { Type_names[nSamples] = RawMatrix[1,nSamples+1];
nSamples = nSamples +1;

}

#Type = vector()
#for (s in 1:NumSamples) { if(Type_names[s] == 'E3'|Type_names[s] =='E2'|Type_names[s] =='E1') {Type[s] = 'E'}}
#for (s in 1:NumSamples) { if(Type_names[s] == 'C3'|Type_names[s] =='C2'|Type_names[s] =='C1') {Type[s] = 'C'}}
#for (s in 1:NumSamples) { if(Type_names[s] == 'C0') {Type[s] = 'C0'}}%%
Batch = vector()
for (s in 1:NumSamples) {Batch[s] = '1'}

Replicates = vector()
NumQC = 0
for (s in 1:NumSamples) {if (Type_names[s]=='QC') {NumQC = NumQC +1}}
NumE4 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E4') {NumE4 = NumE4 +1}}
NumE3 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E3') {NumE3 = NumE3 +1}}
NumE2 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E2') {NumE2 = NumE2 +1}}
NumE1 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E1') {NumE1 = NumE1 +1}}
NumC4 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='C4') {NumC4 = NumC4 +1}}
NumC3 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='C3') {NumC3 = NumC3 +1}}
NumC2 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='C2') {NumC2 = NumC2 +1}}
NumC1 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='C1') {NumC1 = NumC1 +1}}
NumC0 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='C0') {NumC0 = NumC0 +1}}




violin_data = read.xlsx("sigpeaks_EC.xlsx",sheet = 1,colNames= FALSE,rowNames = TRUE)
#lastcol = length(violin_data)
#seclastcol = lastcol - 1

#violin_data = violin_data[,-seclastcol]
#violin_data = violin_data[,-seclastcol]
colnames(violin_data) = Type_names

data_test_violin = t(violin_data[1,] )#take the target data
T_sampletype = t(Type_names)

real_violin = as.data.frame(cbind(Type_names,data_test_violin))
colnames(real_violin) = c("Type","Value")
real_violin$Value = as.numeric(real_violin$Value)
#real_violin$Type = factor(real_violin$Type)
#ggplot(real_violin)+
  #geom_violin(aes(Type,Value,fill=Type,color = Type),width = 1.5)+
  #geom_boxplot(aes(Type,Value),fill = "#a6a7ac",color = "#a6a7ac",width = 0.1,outlier.shape = NA)+
  #scale_fill_manual(values = c("#d1d2d2" "#fbd3b9" "#DADAEB" "#BCBDDC" "#9E9AC8" "#807DBA" "#6A51A3" "#54278F" "#3F007D"))+
  #color_fill_manual(values = c("#d1d2d2" "#fbd3b9" "#DADAEB" "#BCBDDC" "#9E9AC8" "#807DBA" "#6A51A3" "#54278F" "#3F007D"))
  #ylim(9,14)+
  #theme_bw()+
  #theme(panel.grid = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1))+
  #xlab("Sample Type")+
  #ylab("gloged value")

E_begin = NumQC+1
E_end = NumQC+NumE4+NumE3+NumE2+NumE1
C_begin = NumQC+NumE4+NumE3+NumE2+NumE1+1
C_end = NumQC+NumE4+NumE3+NumE2+NumE1+NumC4+NumC3+NumC2+NumC1
C0_begin = NumQC+NumE4+NumE3+NumE2+NumE1+NumC4+NumC3+NumC2+NumC1+1
C0_end = NumSamples
            
E_violin = rbind(real_violin[E_begin:E_end,],real_violin[C0_begin:C0_end,])          
C_violin = real_violin[C_begin:C0_end,]

for ( s in 1: nrow(E_violin)) {if (E_violin[s,1] == 'E4') {E_violin[s,1] = 'T4'}}
for ( s in 1: nrow(E_violin)) {if (E_violin[s,1] == 'E3') {E_violin[s,1] = 'T3'}}
for ( s in 1: nrow(E_violin)) {if (E_violin[s,1] == 'E2') {E_violin[s,1] = 'T2'}}
for ( s in 1: nrow(E_violin)) {if (E_violin[s,1] == 'E1') {E_violin[s,1] = 'T1'}}
for ( s in 1: nrow(E_violin)) {if (E_violin[s,1] == 'C0') {E_violin[s,1] = 'T0'}}


for ( s in 1: nrow(C_violin)) {if (C_violin[s,1] == 'C4') {C_violin[s,1] = 'T4'}}
for ( s in 1: nrow(C_violin)) {if (C_violin[s,1] == 'C3') {C_violin[s,1] = 'T3'}}
for ( s in 1: nrow(C_violin)) {if (C_violin[s,1] == 'C2') {C_violin[s,1] = 'T2'}}
for ( s in 1: nrow(C_violin)) {if (C_violin[s,1] == 'C1') {C_violin[s,1] = 'T1'}}
for ( s in 1: nrow(C_violin)) {if (C_violin[s,1] == 'C0') {C_violin[s,1] = 'T0'}}



exposed_violin = ggplot(E_violin)+
geom_violin(aes(Type,Value,fill=Type,color = Type),width = 1.5)+
geom_boxplot(aes(Type,Value),fill = "#a6a7ac",color = "#a6a7ac",width = 0.1,outlier.shape = NA)+
scale_fill_manual(values = c("#d1d2d2", "#fbd3b9", "#a1c9e5", "#417bb9", "#9E9AC8"))+
scale_color_manual(values = c("#d1d2d2", "#fbd3b9", "#a1c9e5", "#417bb9", "#9E9AC8"))+
  #scale_y_continuous(breaks = seq(8,10,12))+
  ylim(9.5,11.5)+
theme_bw()+
theme(panel.grid = element_blank(),
axis.text.x = element_text(angle = 45,hjust = 1))+
xlab("Time Point")+
ylab("gloged value")+
  annotate("text",x = 1.3, y = 9.95,label = NumC0,color ="#a6a7ac" )+
annotate("text",x = 2.3, y = 10.05,label = NumE1,color ="#fbd3b9" )+
  annotate("text",x = 3.35, y = 10.3,label = NumE2,color ="#a1c9e5" )+
  annotate("text",x = 4.2, y = 10,label = NumE3,color ="#417bb9" )+
  annotate("text",x = 5.45, y = 10.45,label = NumE4,color ="#9E9AC8" )

ggsave(exposed_violin,filename = "exposed_violin.jpeg",width = 16, height = 12,device = NULL)


control_violin = ggplot(C_violin)+
  geom_violin(aes(Type,Value,fill=Type,color = Type),width = 1.5)+
  geom_boxplot(aes(Type,Value),fill = "#a6a7ac",color = "#a6a7ac",width = 0.1,outlier.shape = NA)+
  scale_fill_manual(values = c("#d1d2d2", "#fbd3b9", "#a1c9e5", "#417bb9", "#9E9AC8"))+
  scale_color_manual(values = c("#d1d2d2", "#fbd3b9", "#a1c9e5", "#417bb9", "#9E9AC8"))+
  #scale_y_continuous(breaks = seq(8,10,12))+
  ylim(9.5,11.5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))+
  xlab("Time Point")+
  ylab("gloged value")+
  annotate("text",x = 1.3, y = 9.95,label = NumC0,color ="#a6a7ac" )+
  annotate("text",x = 2.3, y = 9.9,label = NumE1,color ="#fbd3b9" )+
  annotate("text",x = 3.35, y = 9.95,label = NumE2,color ="#a1c9e5" )+
  annotate("text",x = 4.2, y = 9.95,label = NumE3,color ="#417bb9" )+
  annotate("text",x = 5.45, y = 10.05,label = NumE4,color ="#9E9AC8" )

ggsave(control_violin,filename = "control_violin.jpeg",width = 16, height = 12,device = NULL)


