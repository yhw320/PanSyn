library(pheatmap)
## Set working directory
setwd("/PATH/workdirectory")
## Load the matrix file
data <- read.table("Syntenic_percentages.matrix",header=T,row.names=1,sep="\t")
## Open a PDF device and set the filename
pdf("heatmap.pdf")
## Plot the heatmap
pheatmap(data, 
color = colorRampPalette(rev(c("#E35632", "#EDF8C0","#2E79B2")))(100), 
cluster_rows = F,
cluster_cols = F,
fontsize_col = 9,
fontsize_row = 9, 
scale = "none",
border_color = "#c3c3c3",
border_width = 0.05,
breaks = seq(0,1,length.out=101))
## Close the PDF device and save the file
dev.off()
