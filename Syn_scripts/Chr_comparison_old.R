arg <- commandArgs(T)
library(vctrs)
library(RColorBrewer)
library(RIdeogram)
setwd(arg[3])
x <- read.table(arg[2], sep = "\t", header = T, stringsAsFactors = F)
karyotype <-as.data.frame(x)
y <- read.table(arg[1], sep = "\t", header = T, stringsAsFactors = F)
density<-as.data.frame(y)

ideogram(karyotype = karyotype)

col1<-brewer.pal(12,'Paired')
col<-colorRampPalette(c(col1))(n = arg[5])
ideogram(karyotype = karyotype, overlaid = density,colorset1 =col,output="chromosome.svg")
na<-paste(arg[4],"chromosome",sep="_")
convertSVG("chromosome.svg", file=na, device = "pdf",width=12,height=8)
