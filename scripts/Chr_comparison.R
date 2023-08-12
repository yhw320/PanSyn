arg <- commandArgs(T)
library(vctrs)
library(ggplot2)
x<-read.table(arg[2],header=T,sep="\t")
karyotype<-as.data.frame(x)
chrMaxLen <- max(karyotype$End)+100

x1<-read.table(arg[1],header=T,sep="\t")
mergedDF<-as.data.frame(x1)

num1 <-as.numeric(arg[6])
num2 <-as.numeric(arg[7])

chr_data <- read.table(arg[8], header=F,sep="\t")
chr_data2<-as.data.frame(chr_data)

colors <- chr_data2[,2]
myClr <- colors


      ## magic number to draw rectangle
offset = ifelse (max(karyotype$Chr)<=7, 0.05, ifelse (max(karyotype$Chr)<=11, 0.07, ifelse (max(karyotype$Chr)<=16, 0.11, ifelse (max(karyotype$Chr)<=25,0.18,0.3)))) 

	p<-ggplot() +
        geom_segment(data = karyotype,
                     aes(y = Chr, yend = Chr, x = 0, xend = End),
                     lineend = "round", color = "#E6E6E6", linewidth = 6) +
        scale_x_continuous("Length (bp)", n.breaks = 5, limits=c(0,chrMaxLen),labels = scales::comma)+
        scale_y_continuous("Chromosome", breaks = karyotype$Chr, labels = karyotype$Chr ) +
        
        geom_rect(data=mergedDF, mapping=aes(xmin=Start, xmax=End, ymin=as.integer(Chr)-offset, 
                                             ymax=as.integer(Chr)+offset), 
                  fill=myClr[mergedDF$Ancestor_chr], show.legend = FALSE)+
        ##### Show ancestral chromosome number 
        #geom_text(data=mergedDF, aes(x=start+(end-start)/2, y=chr, label=ancestralChr), size=2) +
        coord_flip()  +
        theme(plot.title = element_text(size=15, face="bold"), 
              text = element_text(size=15),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = .5, linetype = "solid")) # enable axis lines #, axis.ticks = element_blank() )

na<-paste(arg[4],"chromosome.pdf",sep="_")


output_dir <- arg[3]
output_file <- file.path(output_dir, na)

# 将图形保存到指定的目录
ggsave(p, file = output_file, width = num1, height = num2)