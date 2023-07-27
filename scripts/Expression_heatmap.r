library(pheatmap)

arg <- commandArgs(T)


df = read.delim(arg[2], 
                header = T,
                sep = "\t",
                row.names = 1,
                fill=T)

wi<-as.numeric(arg[3])
he<-as.numeric(arg[4])
# 设置保存文件的目录
directory <- arg[1]

# 设置文件名
filename <- "Expression_heatmap.pdf"

# 使用 file.path() 函数构建完整的文件路径
file_path <- file.path(directory, filename)

pdf(file=file_path, height=he,width=wi)
#df<-log10(df+1)
pheatmap(df, 
         show_colnames = TRUE,
         show_rownames=TRUE, 
         color = colorRampPalette(c('#448FD7','#ffffff','#FF8A60'))(50), 
         annotation_legend=TRUE,
         border_color="white", 
         scale="column", 
         cluster_rows = FALSE,
         cluster_cols = FALSE
)

dev.off()
