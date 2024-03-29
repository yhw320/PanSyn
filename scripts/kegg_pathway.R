arg <- commandArgs(T)

library("vctrs")
library("ggplot2")

pathway1 = read.table(arg[2],header=T,sep="\t")
p1 = ggplot(pathway1,aes(RichFactor,Description))
p1=p1 + geom_point()
p1=p1 + geom_point(aes(size=Count))
pbubble1 = p1 + geom_point(aes(size=Count,color=-1*LogP))
pbubble1 =pbubble1+ scale_colour_gradient(low="blue",high="red")
pr1 = pbubble1  + 
  labs(color=expression(-log[10](Pvalue)),size="Count",
       x="%",y="Pathway Name",title="")
pr1=pr1 + theme_bw() +theme(axis.text.y =element_text(size=12,color="black"),axis.title.x =element_text(size = 8), title =element_text(size = 10))

wi<-as.numeric(arg[3])
he<-as.numeric(arg[4])


output_dir <- arg[1]
output_file <- file.path(output_dir, "kegg.pdf")

# 将图形保存到指定的目录
ggsave(pr1, file = output_file, width=wi, height=he)