## Alignment with lastz
setwd("/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/try/1")
library(CNEr)


## Config
axtDir <- "/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/try/1"

assemblyQuery <- "/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/inputDir/hg38.2bit"
assemblyTarget <- "/mnt/inspurfs/home/liyl/yuhw/macrosyn-ceshi-data12/Microsyn/supplement/CNEr/inputDir/mm10.2bit"

## Pipeline
lavs <- lastz(assemblyTarget, assemblyQuery, 
              outputDir=axtDir,
	      chrsQuery=c("chr14"),
	      chrsTarget=c("chr12"),
              distance="medium", mc.cores=24)
lavs <- list.files(path=".", pattern="\\.lav$")

psls <- lavToPsl(lavs, removeLav=FALSE, binary="lavToPsl")

chains <- axtChain(psls, assemblyTarget=assemblyTarget,
                   assemblyQuery=assemblyQuery, distance="medium",
                   removePsl=FALSE, binary="axtChain")

allChain <- chainMergeSort(chains, assemblyTarget, assemblyQuery,
              allChain=file.path(axtDir,
                         paste0(sub("\\.2bit$", "", basename(assemblyTarget),
                                    ignore.case=TRUE), ".", 
                                sub("\\.2bit$", "", basename(assemblyQuery), 
                                    ignore.case=TRUE), ".all.chain")),
                           removeChains=FALSE, binary="chainMergeSort")

allPreChain <- chainPreNet(allChain, assemblyTarget, assemblyQuery,
                           allPreChain=file.path(axtDir,
                                      paste0(sub("\\.2bit$", "", 
                                                 basename(assemblyTarget),
                                                 ignore.case = TRUE), ".", 
                                             sub("\\.2bit$", "",
                                                 basename(assemblyQuery),
                                                 ignore.case = TRUE),
                                                 ".all.pre.chain")),
                           removeAllChain=FALSE, binary="chainPreNet")

## Keep the best chain and add synteny information
netSyntenicFile <- chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery,
                     netSyntenicFile=file.path(axtDir,
                                               paste0(sub("\\.2bit$", "",
                                                      basename(assemblyTarget),
                                                      ignore.case = TRUE), ".",
                                                      sub("\\.2bit$", "",
                                                      basename(assemblyQuery),
                                                      ignore.case = TRUE),
                                              ".noClass.net")),
                     binaryChainNet="chainNet", binaryNetSyntenic="netSyntenic")

netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery,
         axtFile=file.path(axtDir,
                           paste0(sub("\\.2bit$", "",
                                      basename(assemblyTarget),
                                      ignore.case = TRUE), ".",
                                  sub("\\.2bit$", "",
                                      basename(assemblyQuery),
                                      ignore.case = TRUE),
                                  ".net.axt")),
             removeFiles=FALSE,
             binaryNetToAxt="netToAxt", binaryAxtSort="axtSort")

