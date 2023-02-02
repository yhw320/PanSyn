## Compile the annotation from refGene and ensGene

setwd("/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/gtf")
library(rtracklayer)

#####mm10

assembly <- "galGal6"
ensGene <- import("galGal6.ensGene.gtf", format="gtf")
refGene <- import("galGal6.refGene.gtf", format="gtf")

filters <- reduce(c(ensGene[ensGene$type=="CDS"], 
                    refGene[refGene$type=="CDS"]), ignore.strand=TRUE)
export.bed(filters, con=paste0("filter_regions.", assembly, ".bed"))

#####hg38

assembly <- "hg38"
ensGene <- import("hg38.ensGene.gtf", format="gtf")
refGene <- import("hg38.refGene.gtf", format="gtf")

filters <- reduce(c(ensGene[ensGene$type=="CDS"], 
                    refGene[refGene$type=="CDS"]), ignore.strand=TRUE)
export.bed(filters, con=paste0("filter_regions.", assembly, ".bed"))
