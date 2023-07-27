## makeGRBs for sear urchin

args <- commandArgs(TRUE)


library(CNEr)
library(GenomicRanges)
library(rtracklayer)
library(RSQLite)
load(args[2])

# 保存当前的临时文件夹路径
old_tempdir <- tempdir()

# 指定新的临时文件夹路径
new_tempdir <- arg[1]

# 在特定代码段中使用新的临时文件夹路径
tempdir(new_tempdir)
# 在此处编写需要使用新临时文件夹的代码


dbName <- "cne.sqlite"
readAncoraIntoSQLite(args[3],
                     dbName=dbName,overwrite=TRUE)

mydb <- dbConnect(RSQLite::SQLite(), "cne.sqlite")

ar1<-as.numeric(args[4])
ar2<-as.numeric(args[5])


query <- sprintf("SELECT * FROM SpeA_SpeB_%s_%s", ar1, ar2)
cnes <- dbGetQuery(mydb, query)

#cnes <- dbGetQuery(mydb, 'SELECT * FROM SpeA_SpeB_48_50')
cnesGR<- GRanges(seqnames=cnes$first.seqnames,
                  ranges=IRanges(start=cnes$first.start,
                                 end=cnes$first.end),
                  strand="*",
                  seqinfo=seqinfo(TwoBitFile(args[6]))
                  )
grangesList <- GRangesList(Species=cnesGR)
ratio1<-as.numeric(args[7])
ratio <- ratio1

win<-as.numeric(args[8])
min<-as.numeric(args[9])



outputDir <- args[1]

grbsSpecies <- makeGRBs(grangesList, winSize=win, genes=NULL, ratio=ratio,
                         background="genome", minCNEs=min)

export.bed(grbsSpecies, paste0(outputDir, "grbsSpeA-SpeB_", ar1 ,"_", ar2, ".bed"))

# 恢复默认的临时文件夹路径
tempdir(old_tempdir)
