

arg <- commandArgs(T)
library(vctrs)
library(ggplot2)
setwd(arg[1])

x<-read.table(arg[2],header=T,sep="\t")
data<-as.data.frame(x)

p<-ggplot(data,mapping=aes(x=Categories,y=Values,fill=Categories))+geom_bar(stat="identity",position=position_dodge(0.6),width=0.6) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + scale_y_continuous(expand=c(0,0),breaks = NULL) + theme(axis.text.x=element_text(size=8,hjust=0.5,vjust=0.5,color= "black"),axis.title.x=element_text(size=12,color= "black"),axis.text.y=element_text(size=10,color= "black"),axis.title.y=element_text(size=12,color= "black"))+xlab(NULL)+ylab(NULL)+theme(legend.position="none")

ggsave(p, file="Categories.pdf", width=10, height=4)
