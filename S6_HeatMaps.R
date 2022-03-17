#BiocManager::install("ComplexHeatmap")
#BiocManager::install("circlize")
#BiocManager::install("grid")
suppressPackageStartupMessages({library (ggplot2);
library(openxlsx);
library(ggthemes)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
  library(grid)
})
setwd("C:/Users/10938/Desktop/Jacobson Lab/2021.10-/4 time points")

data <- read.xlsx("Zscored_sigdata.xlsx",sheet = 1,colNames= TRUE,rowNames = FALSE)
data <- as.matrix(data)
names <- c("QCs","E4","E3","E2","E1","C4","C3","C2","C1","C0")
data <- t(data)
heatmap(data)
names = t(names)
colnames(data) = names
col_fun <- colorRamp2(c(-2.5,0,2.5),c("#ff0000","#ff8000","#ffff00"))

Heatmap(data,
        #color
        col = col_fun,
        #height of the dendrogram
        rect_gp = gpar(col = "white", lwd = 1),
        
        row_dend_width = unit(2,"cm"),
        #size of labels
        row_names_gp = gpar(fontsize = 7 , fontface = "italic"),
        column_names_gp = gpar(fontsize= 7 ),
        #note in the cells
        #cell_fun = function(j,i,x,y,width,height,fill) {
          #grid.text(sprintf("%.1f",data[i,j]),x,y,gp = gpar(fontsize = 5))},
        )
