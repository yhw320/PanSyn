arg <- commandArgs(T)

library("vctrs")
library("ggplot2")

data<-read.table(arg[3],header=T,sep="\t")
data3<-as.data.frame(data)

x<-read.table(arg[2],header=F,sep="\t")
x<-as.data.frame(x)

p1<- ggplot(data3,aes(Chr,Number,fill=Ancestor))+
     geom_bar(stat="identity",position="stack",width=0.5)+
     scale_fill_manual(breaks=c(x$V1),values=c(x$V2))+
     theme_bw()+
     theme(axis.ticks.length=unit(0.5,'cm'))
p2 <- p1+ theme(axis.line = element_blank(), axis.text.y=element_blank(),axis.text.x=element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),panel.border=element_blank(),panel.background = element_rect(fill = NA),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
   labs(title = "Chromosome Rearrangements",x = NULL, y = NULL)+
   theme(axis.text.x = element_text(size = 20),plot.title = element_text(size = 25))+
   theme(legend.title = element_text(size = 15,face = "bold"))+
   theme(legend.key.size = unit(0.2, "inches"))
   #+geom_hline(yintercept = 1,linetype=5,col="red")   
na<-paste(arg[4],"chromsome_break_fusion.pdf",sep="-")


output_dir <- arg[1]
output_file <- file.path(output_dir, na)

# 将图形保存到指定的目录
ggsave(p2, file = output_file, width = 22, height = 4)