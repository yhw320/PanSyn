# Make horizonplot
library(CNEr)
library(Gviz)
library(GenomicFeatures)
library(biomaRt)
library(rtracklayer)


setwd("C:/Users/HP/Desktop")
load("work.RData")
## prepare CNE.sqlite
dbName <- "cne.sqlite"
readAncoraIntoSQLite("cne2wBf_Hg38_Mm10_35_50",
                     dbName=dbName)
readAncoraIntoSQLite("cne2wBf_Hg38_Mm10_40_50",
                     dbName=dbName)
readAncoraIntoSQLite("cne2wBf_Hg38_Mm10_45_50",
                     dbName=dbName)


## Gene tracks for Gviz
genome <- "hg38"
chr <- "chr14"
start <- 35080000
end <- 39080000
axisTrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome=genome, chr=chr)

## CNE density tracks for Gviz
windowSize <- 300L
minLength <- 50L
cne1 <- CNEDensity(dbName=dbName,tableName="Hg38_Mm10_35_50",
             whichAssembly="first", chr=chr, start=start,
             end=end, windowSize=windowSize,
             minLength=minLength)

cne2 <- CNEDensity(dbName=dbName, tableName="Hg38_Mm10_40_50",
                   whichAssembly="first", chr=chr, start=start,
                   end=end, windowSize=windowSize,
                   minLength=minLength)

cne3 <- CNEDensity(dbName=dbName, tableName="Hg38_Mm10_45_50",
                   whichAssembly="first", chr=chr, start=start,
                   end=end, windowSize=windowSize,
                   minLength=minLength)

horizon.scale <- max(cne1$score)/3

## CNE tracks for Gviz
dTrack1 <- DataTrack(range=cne1,
                     genome=genome, type="horiz",
                     horizon.scale=horizon.scale,
                     fill.horizon=c("#B41414", "#E03231", "#F7A99C",
                                    "yellow", "orange", "red"),
                     name="MM10\n35/50", background.title="brown")

dTrack2 <- DataTrack(range=cne2,
                     genome=genome, type="horiz",
                     horizon.scale=horizon.scale,
                     fill.horizon=c("#B41414", "#E03231", "#F7A99C",
                                    "yellow", "orange", "red"),
                     name="MM10\n40/50", background.title="brown")

dTrack3 <- DataTrack(range=cne3,
                     genome=genome, type="horiz",
                     horizon.scale=horizon.scale,
                     fill.horizon=c("#B41414", "#E03231", "#F7A99C",
                                    "yellow", "orange", "red"),
                     name="MM10\n45/50", background.title="brown")
options(ucscChromosomeNames=FALSE)
CNEr:::savefig("Hg38-Mm10", colormodel="cmyk", height=5, width=15,
               type="pdf", family="sans")
plotTracks(list(axisTrack, 
                dTrack1, dTrack3), from=start, to=end)
dev.off()

save.image("work-Hg38-Mm10.RData")


