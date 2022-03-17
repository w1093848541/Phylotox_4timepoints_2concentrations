# install BiocManager if not present
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# install structToolbox and dependencies
#BiocManager::install("structToolbox")
#BiocManager::install("mixOmics")
#a
#BiocManager::install(c('pmp', 'ropls', 'BiocFileCache'))
#install.packages(c('cowplot', 'openxlsx'))
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
MATLABresult = read.xlsx("C:\\Users\\10938\\Desktop\\Jacobson Lab\\2021.10-\\4 time points\\Combined_withQC_4tp_all.xlsx",sheet = 1,colNames= FALSE,rowNames = FALSE);
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
noQC = read.xlsx("C:\\Users\\10938\\Desktop\\Jacobson Lab\\2021.10-\\4 time points\\Combined_withQC_4tp_all.xlsx",sheet = 1,colNames= TRUE,rowNames = TRUE);
noQC = as.matrix(noQC)
colnames(noQC) = Type_names
rownames(noQC) = Mass_list
Type = Type_names
OPLS_DA_GROUPS = factor(Type_names)
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
for (s in 1:NumQC) {Replicates[s] = as.character(s)}
for (s in 1:NumE4) {Replicates[s+NumQC] = as.character(s)}
for (s in 1:NumE3) {Replicates[s+NumE4+NumQC] = as.character(s)}
for (s in 1:NumE2) {Replicates[s+NumE4+NumE3+NumQC] = as.character(s)}
for (s in 1:NumE1) {Replicates[s+NumE4+NumE3+NumE2+NumQC] = as.character(s)}
for (s in 1:NumC4) {Replicates[s+NumE4+NumE3+NumE2+NumE1+NumQC] = as.character(s)}
for (s in 1:NumC3) {Replicates[s+NumE4+NumE3+NumE2+NumE1+NumC4+NumQC] = as.character(s)}
for (s in 1:NumC2) {Replicates[s+NumE4+NumE3+NumE2+NumE1+NumC4+NumC3+NumQC] = as.character(s)}
for (s in 1:NumC1) {Replicates[s+NumE4+NumE3+NumE2+NumE1+NumC4+NumC3+NumC2+NumQC] = as.character(s)}
for (s in 1:NumC0) {Replicates[s+NumE4+NumE3+NumE2+NumE1+NumC4+NumC3+NumC2+NumC1+NumQC] = as.character(s)}
PhylotoxData@colData@listData[["Sample_Rep"]] = Replicates
PhylotoxData@colData@listData[["Batch"]] = Batch
PhylotoxData@colData@listData[["Class"]] = Type
PhylotoxData@colData@nrows = as.integer(NumSamples)
PhylotoxData@elementMetadata@nrows = as.integer(NumMetabolites)
PhylotoxData@colData@listData[["Class2"]] = Type_names
PhylotoxData@assays@data@listData[[1]] = noQC
PhylotoxData@NAMES= Mass_list
PhylotoxData@colData@rownames = Type_names
DE = as.DatasetExperiment(PhylotoxData)
DE$sample_meta$run_order = 1:nrow(DE)
Type=as.character(DE$sample_meta$Class)
Type[Type != 'QC'] = 'Sample'
DE$sample_meta$Type = factor(Type)
DE$sample_meta$Batch = factor(DE$sample_meta$Batch)
DE$sample_meta$Class = factor(DE$sample_meta$Class)





M5 = pqn_norm(qc_label='QC',factor_name='Type') 
  
M5 = model_apply(M5,DE)

PQN_result = M5@normalised@value@assays@data@listData[[1]]
PQN_result = t(PQN_result)
colnames(PQN_result) = Type_names
write.csv(PQN_result,file = "C:\\Users\\10938\\Desktop\\Jacobson Lab\\2021.10-\\4 time points\\PQN_predictedresult.csv")

M6 = knn_impute(neighbours=5,by='samples') +
  glog_transform(qc_label='QC',factor_name='Type')

M6 = model_apply(M6,predicted(M5))

GlogedData = predicted(M6)
glog_result = GlogedData@assays@data@listData[[1]]
glog_result = t(glog_result)
colnames(glog_result) = Type_names
write.csv(glog_result,file = "C:\\Users\\10938\\Desktop\\Jacobson Lab\\2021.10-\\4 time points\\glog_predictedresult.csv")

# PCA
M7  = mean_centre() + PCA(number_components = 2)

# apply model sequence to data
M7 = model_apply(M7,predicted(M6))

# plot pca scores
C = pca_scores_plot(factor_name=c('Sample_Rep','Class'),ellipse='none')
chart_plot(C,M7[2]) + coord_fixed() +guides(colour=FALSE) + scale_shape_manual(values = 15:24)

PhylotoxPCA = M7@models[[2]]@scores@value@assays@data@listData[[1]]
write.csv(PhylotoxPCA,file = "C:\\Users\\10938\\Desktop\\Jacobson Lab\\2021.10-\\4 time points\\WithQCs_pca.csv")

PhylotoxData2 = MTBLS79

nRemove = 1
while(nRemove <= NumQC){Replicates = Replicates[-1];
Batch = Batch[-1];
Type = Type[-1];
Type_names = Type_names[-1];
M6@models[[2]]@transformed@value@assays@data@listData[[1]]= M6@models[[2]]@transformed@value@assays@data@listData[[1]][-1,];
nRemove = nRemove + 1;
}

PhylotoxData2@colData@listData[["Sample_Rep"]] = Replicates
PhylotoxData2@colData@listData[["Batch"]] = Batch
PhylotoxData2@colData@listData[["Class"]] = Type
PhylotoxData2@colData@nrows = as.integer(NumSamples-NumQC)
PhylotoxData2@elementMetadata@nrows = as.integer(NumMetabolites)
PhylotoxData2@colData@listData[["Class2"]] = Type_names
PhylotoxData2@assays@data@listData[[1]] = t(M6@models[[2]]@transformed@value@assays@data@listData[[1]])
PhylotoxData2@NAMES= Mass_list
PhylotoxData2@colData@rownames = Type_names
DE = as.DatasetExperiment(PhylotoxData2)
DE$sample_meta$run_order = 1:nrow(DE)
Type=as.character(DE$sample_meta$Class)
Type[Type != 'QC'] = 'Sample'
DE$sample_meta$Type = factor(Type)
DE$sample_meta$Batch = factor(DE$sample_meta$Batch)
DE$sample_meta$Class = factor(DE$sample_meta$Class)



# PCA
M7  = mean_centre() + PCA(number_components = 2)

# apply model sequence to data
M7 = model_apply(M7,DE)

# plot pca scores
C = pca_scores_plot(factor_name=c('Sample_Rep','Class2'),ellipse='none')
chart_plot(C,M7[2]) + coord_fixed() +guides(colour=FALSE) + scale_shape_manual(values = 16:24)

PhylotoxPCA = M7@models[[2]]@scores@value@assays@data@listData[[1]]
write.csv(PhylotoxPCA,file = "C:\\Users\\10938\\Desktop\\Jacobson Lab\\2021.10-\\3concentrations\\200_ppmf\\WithQCsbutremoved_pca.csv")

gloged = read.csv("glog_predictedresult.csv",header = TRUE,row.names = 1)

gloged = t(gloged)

plsda.datatm <- plsda(gloged,OPLS_DA_GROUPS,ncomp = 3)
jpeg(filename = "plada_result.jpeg",width = 5, height = 5, units = "in",res =300 )

plotIndiv(plsda.datatm,ind.names = FALSE,legend = TRUE, ellipse = FALSE,title = "PLS-DA - Result")
dev.off()


background <- background.predict(plsda.datatm, comp.predicted = 2, dist = "max.dist")
jpeg(filename = "plsda_max.jpeg",width = 5, height = 5, units = "in",res =300)

 plotIndiv(plsda.datatm, comp = 1:2, ind.names = FALSE, title = "Maximum Distance",
          legend = TRUE, background = background, ellipse = FALSE)
 dev.off()
#opls(gloged,OPLS_DA_GROUPS,predI = 1)