arg <- commandArgs(T)
setwd(arg[4])

library(vctrs)
library(ggstance)
library(ggplot2)
library(ggtree)


tree<-read.tree(arg[3])
data=fortify(tree)

x<-read.table(arg[1],header=T,sep="\t")
mydata<-as.data.frame(x)

x<-read.table(arg[2],header=T,sep="\t")
mydata2<-as.data.frame(x)

group_file <- read.table(arg[5],header = T,row.names = 1)

groupInfo <- split(row.names(group_file), group_file$Group)
tree <- groupOTU(tree, groupInfo)

p <- ggtree(tree)+
geom_tiplab(size=7,color="black") +
theme_tree2()+
xlim_expand(c(0,max(data$x)*1.2),'Tree')



#ggThemeAssistGadget(p2)

p2 <- facet_plot(p, panel="Barplot", data=mydata, geom = geom_segment,
                 mapping = aes(x =0,xend= Total_ancient_gene_families,y=y+0.1,yend=y+0.1),size=4,color="#20BAB2")+
		 xlim_expand(c(0,1),'Barplot')


p3 <- facet_plot(p2, panel="Barplot", data=mydata, geom = geom_segment,
                 mapping = aes(x =0,xend= Total_gene_families,y=y-0.1,yend=y-0.1),size=4,color="#FF728A")+
		 xlim_expand(c(0,1),'Barplot')

p4 <- facet_plot(p3, panel="Stacked Barplot", data=mydata2, geom = geom_barh,col="white",
                 mapping = aes(x = Number, fill=Range), 
                 stat='identity',width=.6)+scale_fill_manual(values=c("#F180B8","#88CBFF","#FFA874","#FACCCE"))

p5 <- facet_plot(p4,panel="Stacked Barplot", data=mydata2,geom= geom_text,col="grey20",mapping = aes(label=round(Number,digits=0),x=round(Position,digits=0)))

p6 <- facet_widths(p5,widths=c(1.2,1,1.5))

num1 <-as.numeric(arg[6])
num2 <-as.numeric(arg[7])

ggsave(p6, file="Ancient_gene_families.pdf",width=num1,height=num2)

