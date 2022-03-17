
library (ggplot2);
library(openxlsx);
library("ggthemes")
library("RColorBrewer")
setwd("C:/Users/10938/Desktop/Jacobson Lab/2021.10-/4 time points")
T1data = read.xlsx("C:\\Users\\10938\\Desktop\\Jacobson Lab\\2021.10-\\4 time points\\qvalues_allT1.xlsx",sheet = 1,colNames= FALSE,rowNames = FALSE);



data = T1data[,3:4];
colnames(data) = c("q","log2FoldChange")


data$label <- T1data[,1]
n = nrow(data)
for (q in 1:n) {if (data[q,1] > 0.05) {data[q,3] = NaN}}


q1 <- ggplot(data,aes(log2FoldChange,-log10(q))) +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  geom_vline(xintercept = c(-1.2,1.2),linetype = "dashed")+
  geom_point(aes(size = -log10(q),color = -log10(q))) +
  theme_bw() +
theme(panel.grid = element_blank(),
      legend.position = c(0.99,0.99),
      legend.justification = c(0.99,0.99))  +
  scale_colour_gradientn(values = seq(0,1,0.2), colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  guides(col = guide_colourbar(title = "-log10 q-value"),
         size = "none")+
         #size = guide_legend(title = "-log10 q-value"))
 xlab("log2FoldChange(Exposed/Control)")+
   ylab ("-log10(FDR q-value)")+
  geom_text(aes(label = label,color = -log10(q)),size = 3, vjust = 1.5, hjust = 1)




ggsave(q1,filename = "T1_volcano.png",width = 16, height = 12)


T2data = read.xlsx("qvalues_allT2.xlsx",sheet = 1,colNames= FALSE,rowNames = FALSE);



data = T2data[,3:4];
colnames(data) = c("q","log2FoldChange")


data$label <- T2data[,1]
n = nrow(data)
for (q in 1:n) {if (data[q,1] > 0.05) {data[q,3] = NaN}}


q2<-ggplot(data,aes(log2FoldChange,-log10(q))) +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  geom_vline(xintercept = c(-1.2,1.2),linetype = "dashed")+
  geom_point(aes(size = -log10(q),color = -log10(q))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.99,0.99),
        legend.justification = c(0.99,0.99))  +
  scale_colour_gradientn(values = seq(0,1,0.2), colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  guides(col = guide_colourbar(title = "-log10 q-value"),
         size = "none")+
  #size = guide_legend(title = "-log10 q-value"))
  xlab("log2FoldChange(Exposed/Control)")+
  ylab ("-log10(FDR q-value)")+
  geom_text(aes(label = label,color = -log10(q)),size = 3, vjust = 1.5, hjust = 1)




ggsave(q2,filename = "T2_volcano.png",width = 16, height = 12)

T3data = read.xlsx("qvalues_allT3.xlsx",sheet = 1,colNames= FALSE,rowNames = FALSE);



data = T3data[,3:4];
colnames(data) = c("q","log2FoldChange")


data$label <- T3data[,1]
n = nrow(data)
for (q in 1:n) {if (data[q,1] > 0.05) {data[q,3] = NaN}}


q3<-ggplot(data,aes(log2FoldChange,-log10(q))) +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  geom_vline(xintercept = c(-1.2,1.2),linetype = "dashed")+
  geom_point(aes(size = -log10(q),color = -log10(q))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.99,0.99),
        legend.justification = c(0.99,0.99))  +
  scale_colour_gradientn(values = seq(0,1,0.2), colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  guides(col = guide_colourbar(title = "-log10 q-value"),
         size = "none")+
  #size = guide_legend(title = "-log10 q-value"))
  xlab("log2FoldChange(Exposed/Control)")+
  ylab ("-log10(FDR q-value)")+
  geom_text(aes(label = label,color = -log10(q)),size = 3, vjust = 1.5, hjust = 1)




ggsave(q3,filename = "T3_volcano.png",width = 16, height = 12)

T4data = read.xlsx("qvalues_allT4.xlsx",sheet = 1,colNames= FALSE,rowNames = FALSE);



data = T4data[,3:4];
colnames(data) = c("q","log2FoldChange")


data$label <- T4data[,1]
n = nrow(data)
for (q in 1:n) {if(data[q,1] > 0.05) {data[q,3] = NaN}}



q4<-ggplot(data,aes(log2FoldChange,-log10(q))) +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  geom_vline(xintercept = c(-1.2,1.2),linetype = "dashed")+
  geom_point(aes(size = -log10(q),color = -log10(q))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.99,0.99),
        legend.justification = c(0.99,0.99))  +
  scale_colour_gradientn(values = seq(0,1,0.2), colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  guides(col = guide_colourbar(title = "-log10 q-value"),
         size = "none")+
  #size = guide_legend(title = "-log10 q-value"))
  xlab("log2FoldChange(Exposed/Control)")+
  ylab ("-log10(FDR q-value)")+
  geom_text(aes(label = label,color = -log10(q)),size = 3, vjust = 1.5, hjust = 1)




ggsave(q4,filename = "T4_volcano.png",width = 16, height = 12)

