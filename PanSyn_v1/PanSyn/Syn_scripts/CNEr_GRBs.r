## makeGRBs for sear urchin
setwd("/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/try/4")

library(CNEr)
library(GenomicRanges)
library(rtracklayer)
library(RSQLite)
load("work-Hg38-Mm10.RData")

dbName <- "cne.sqlite"
readAncoraIntoSQLite("/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/try/3/Ancora/cne2wBf_Hg38_Mm10_35_50",
                     dbName=dbName)

mydb <- dbConnect(RSQLite::SQLite(), "cne.sqlite")
cnes <- dbGetQuery(mydb, 'SELECT * FROM Hg38_Mm10_35_50')
cnesGR<- GRanges(seqnames=cnes$first.seqnames,
                  ranges=IRanges(start=cnes$first.start,
                                 end=cnes$first.end),
                  strand="*",
                  seqinfo=seqinfo(TwoBitFile("/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/inputDir/hg38.2bit"))
                  )
grangesList <- GRangesList(Hg38=cnesGR)
ratio <- 0.95
grbsHg38 <- makeGRBs(grangesList, winSize=300, genes=NULL, ratio=ratio,
                         background="genome", minCNEs=10)
export.bed(grbsHg38, paste0("grbsHg38-", ratio, ".bed"))
