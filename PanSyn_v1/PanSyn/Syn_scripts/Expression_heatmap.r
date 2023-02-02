library(pheatmap)

setwd("C:/Users/HP/Desktop")

df = read.delim("1.txt", 
                header = T,
                sep = "\t",
                row.names = 1,
                fill=T)
pdf(file='test.pdf', height=3,width=6)
#df<-log10(df+1)
pheatmap(df, 
         show_colnames = TRUE,
         show_rownames=TRUE, 
         fontsize=10,
         color = colorRampPalette(c('#448FD7','#ffffff','#FF8A60'))(50), 
         annotation_legend=TRUE,
         border_color="white", 
         scale="column", 
         cluster_rows = FALSE,
         cluster_cols = FALSE
)

dev.off()
