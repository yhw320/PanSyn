## Compile the annotation from refGene and ensGene

library(rtracklayer)

args <- commandArgs(TRUE)


ensGene_file <- args[2]
refGene_file <- args[3]

ensGene <- import(ensGene_file, format="gtf")
refGene <- import(refGene_file, format="gtf")

filters <- reduce(c(ensGene[ensGene$type=="CDS"], 
                    refGene[refGene$type=="CDS"]), ignore.strand=TRUE)


outputDir <- args[1]
export.bed(filters, con = file.path(outputDir, paste0("filter_regions", ".bed")))


