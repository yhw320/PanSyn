arg <- commandArgs(T)

library("vctrs")
library("ggplot2")
data=read.table(arg[2],header=T,sep="\t")
COLS<-c("#90BFFA","#A0FFA0","#F8A0A0")
GO_Term_order=factor(as.integer(rownames(data)),labels=data$Description)
p1 <- ggplot(data=data,aes(x=GO_Term_order,y=-LogP, fill=Category))+
      geom_bar(stat="identity",width=0.5)+
      scale_fill_manual(values = COLS)+
      coord_flip() + 
      theme_bw()+
      labs(title="The Most Enriched GO Terms",y="-log10(Pvalue)")
p2 <- p1+ theme(axis.line = element_line(size = 0), axis.ticks = element_line(size = 0),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_rect(fill = NA))

wi<-as.numeric(arg[3])
he<-as.numeric(arg[4])


output_dir <- arg[1]
output_file <- file.path(output_dir, "go.pdf")

# 将图形保存到指定的目录
ggsave(p2, file = output_file, width=wi, height=he)