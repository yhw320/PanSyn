arg <- commandArgs(T)
library(vctrs)
library(ggplot2)
setwd(arg[3])
x<-read.table(arg[2],header=T,sep="\t")
karyotype<-as.data.frame(x)
chrMaxLen <- max(karyotype$End)+100

x1<-read.table(arg[1],header=T,sep="\t")
mergedDF<-as.data.frame(x1)

num1 <-as.numeric(arg[6])
num2 <-as.numeric(arg[7])

myClr <- c("#A5CDE1", "#4A91C1", "#569DA3", "#AAD383", "#52AE48", "#899C5A", "#F48988","#E63333", "#EE6B46", "#FAB15B", "#F5861F", "#E19B78", "#B093C4", "#70459A", "#C6B598", "#E5CA75", "#B05A28" )  
      
      ## magic number to draw rectangle
      offset = ifelse (max(karyotype$Chr)<=7, 0.05, ifelse (max(karyotype$Chr)<=11, 0.07, ifelse (max(karyotype$Chr)<=16, 0.11, ifelse (max(karyotype$Chr)<=25,0.18,0.3)))) 

	p<-ggplot() +
        geom_segment(data = karyotype,
                     aes(y = Chr, yend = Chr, x = 0, xend = End),
                     lineend = "round", color = "#E6E6E6", size = 6) +
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
              axis.line = element_line(colour = "black", size = .5, linetype = "solid")) # enable axis lines #, axis.ticks = element_blank() )

na<-paste(arg[4],"chromosome.pdf",sep="_")
ggsave(p, file=na, width = num1, height = num2)
