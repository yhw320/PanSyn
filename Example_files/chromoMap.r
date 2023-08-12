library(chromoMap)
## Set working directory
setwd("/PATH/workdirectory")
## Load three files
x1<-read.table("rheMac3.chromosome",header=F,sep="\t")
x1<-as.data.frame(x1)
##
x2<-read.table("rheMac3_annotation_pos.txt",header=F,sep="\t")
x2<-as.data.frame(x2)
##
x3<-read.table("color.txt",header=F,sep="\t")
x3<-as.data.frame(x3) 
## Plot
chromoMap(list(x1),list(x2),segment_annotation=T,data_based_color_map=T,chr_color=c("#F7F7F7"),data_type="categorical",data_colors = list(x3$V2),discrete.domain = list(x3$V1),legend=T,lg_x=0,lg_y = 750,export.options=T)
