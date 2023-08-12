## Make horizonplot 
library(CNEr)
library(Gviz)
library(GenomicFeatures)
library(biomaRt)
library(rtracklayer)

## Set working directory
setwd("PATH/directory")
## prepare CNE.sqlite
dbName <- "cne1.sqlite"
## Read Data
readAncoraIntoSQLite("cne2wBf_SpeA_SpeB_35_50",dbName=dbName)
readAncoraIntoSQLite("cne2wBf_SpeA_SpeB_45_50",dbName=dbName)
## Specify the chromosome and start/end positions to be plotted
genome <- "HSap"
chr <- "HSap_chr14"
start <- 35080000
end <- 39080000
axisTrack <- GenomeAxisTrack()
options(ucscChromosomeNames=FALSE)
## CNE density tracks for Gviz
windowSize <- 300L
minLength <- 50L
## Prepare CNE1
cne1 <- CNEDensity(dbName=dbName, tableName="SpeA_SpeB_35_50",
             whichAssembly="first", chr=chr, start=start,
             end=end, windowSize=windowSize,
             minLength=minLength)
## Prepare CNE2
cne2 <- CNEDensity(dbName=dbName,
             tableName="SpeA_SpeB_45_50",
             whichAssembly="first", chr=chr, start=start,
             end=end, windowSize=windowSize,
             minLength=minLength)
## CNE tracks for Gviz
dTrack1 <- DataTrack(range=cne1,
                     genome=genome, type="horiz",
                     horizon.scale=max(cne1$score)/3,
                     fill.horizon=c("#B41414", "#E03231", "#F7A99C",
                                    "yellow", "orange", "red"),
                     name="SpeB\n35/50", background.title="brown")
dTrack2 <- DataTrack(range=cne2,
                     genome=genome, type="horiz",
                     horizon.scale=max(cne2$score)/3,
                     fill.horizon=c("#B41414", "#E03231", "#F7A99C",
                                    "yellow", "orange", "red"),
                     name="SpeB\n45/50", background.title="brown")
## 
options(ucscChromosomeNames=FALSE)
## Save file
CNEr:::savefig("SpeA-SpeB", colormodel="rgb", height=7, width=15, type="pdf", family="sans")
##
plotTracks(list(axisTrack, dTrack1, dTrack2), collapseTranscripts=TRUE, transcriptAnnotation="symbol", from=start, to=end)
##
dev.off()
##
save.image("work-SpeA-SpeB.RData")
