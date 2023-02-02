## CNEr pipeline
library(CNEr)
library(rtracklayer)
setwd("/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/try2/2mouse")

cneMm10Hg38 <- CNE(
  assembly1Fn="/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/inputDir/mm10.2bit",
  assembly2Fn="/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/inputDir/hg38.2bit",
  axt12Fn="/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/try2/try1/mm10.hg38.net.axt", 
	axt21Fn="/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/try2/try2/hg38.mm10.net.axt",
  cutoffs1=4L, cutoffs2=4L)
Hg38Filters <- import.bed("/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/gtf/filter_regions.hg38.bed")
Mm10Filters <- import.bed("/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/gtf/filter_regions.mm10.bed")

identities <- c( 35, 40, 45)
windows <-    c( 50, 50, 50)
cneListMm10Hg38 <- ceScan(x=cneMm10Hg38, tFilter=Mm10Filters,
                              qFilter=Hg38Filters,
                              window=as.integer(windows), 
                              identity=as.integer(identities))

cneMergedListMm10Hg38 <- lapply(cneListMm10Hg38, cneMerge)
cneFinalListMm10Hg38 <- lapply(cneMergedListMm10Hg38, blatCNE)

## CNE numbers
lengths(lapply(cneFinalListMm10Hg38, CNEFinal))

# Make tracks of bed and bedgraph
library(doParallel)
library(foreach)
registerDoParallel(cores=4)
foreach(threshold=names(cneFinalListMm10Hg38)) %dopar% {
  makeCNEDensity(CNEFinal(cneFinalListMm10Hg38[[threshold]]),
                 genomeFirst="Mm10",
                 genomeSecond="Hg38",
                 threshold=threshold, windowSizeFirst=30, windowSizeSecond=30)
}

## Output data
save.image(file="work.RData")
save(cneFinalListMm10Hg38, file="cneFinalListMm10Hg38.rda")

### Make Ancora format file
for(threshold in names(cneFinalListMm10Hg38)){
  makeAncoraFiles(CNEFinal(cneFinalListMm10Hg38[[threshold]]),
                  outputDir="Ancora", genomeFirst="Mm10",
                  genomeSecond="Hg38", threshold=threshold)
}

