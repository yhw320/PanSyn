arg <- commandArgs(T)

## CNEr pipeline
library(CNEr)
library(rtracklayer)

cut1<-as.integer(arg[6])
cut2<-as.integer(arg[7])
core<-as.numeric(arg[12])
den<-as.numeric(arg[13])

# 保存当前的临时文件夹路径
old_tempdir <- tempdir()

# 指定新的临时文件夹路径
new_tempdir <- arg[1]

# 在特定代码段中使用新的临时文件夹路径
tempdir(new_tempdir)
# 在此处编写需要使用新临时文件夹的代码


outputDir <- arg[1]

cne_SpeA_SpeB <- CNE(
  assembly1Fn=arg[2],
  assembly2Fn=arg[3],
  axt12Fn=arg[4], 
  axt21Fn=arg[5],
  cutoffs1=cut1, cutoffs2=cut2)
SpeA_Filters <- import.bed(arg[8])
SpeB_Filters <- import.bed(arg[9])




identities <- as.integer(arg[10])
windows <- as.integer(arg[11])

cneList_SpeA_SpeB <- ceScan(x=cne_SpeA_SpeB, tFilter=SpeA_Filters, tSizes=NULL, qSizes=NULL,
                              qFilter=SpeB_Filters,
                              window=as.integer(windows), 
                              identity=as.integer(identities))

cneMergedList_SpeA_SpeB <- lapply(cneList_SpeA_SpeB, cneMerge)
cneFinalList_SpeA_SpeB <- lapply(cneMergedList_SpeA_SpeB, blatCNE)

## CNE numbers
lengths(lapply(cneFinalList_SpeA_SpeB, CNEFinal))

# Make tracks of bed and bedgraph
library(doParallel)
library(foreach)
registerDoParallel(cores=core)
foreach(threshold=names(cneFinalList_SpeA_SpeB)) %dopar% {
  makeCNEDensity(CNEFinal(cneFinalList_SpeA_SpeB[[threshold]]),
                 genomeFirst="SpeA",
                 genomeSecond="SpeB",
                 threshold=threshold, windowSizeFirst=den, windowSizeSecond=den)
}
## Output data
z <- paste("work", identities, windows, sep="_")
z1 <- paste(z, ".RData", sep="")


save.image(file = file.path(outputDir, z1))
save(cneFinalList_SpeA_SpeB, file = file.path(outputDir, "cneFinalList_SpeA_SpeB.rda"))

### Make Ancora format file
for (threshold in names(cneFinalList_SpeA_SpeB)) {
  makeAncoraFiles(CNEFinal(cneFinalList_SpeA_SpeB[[threshold]]),
                  outputDir = file.path(outputDir, "Ancora"),
                  genomeFirst = "SpeA",
                  genomeSecond = "SpeB", threshold = threshold)
}

# 恢复默认的临时文件夹路径
tempdir(old_tempdir)
