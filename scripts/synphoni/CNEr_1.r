arg <- commandArgs(T)


## Alignment with lastz
library(CNEr)


## Config
axtDir <- arg[1]



assemblyTarget <- arg[2]
assemblyQuery <- arg[3]
core<-as.numeric(arg[7])

# 保存当前的临时文件夹路径
old_tempdir <- tempdir()

# 指定新的临时文件夹路径
new_tempdir <- arg[1]

# 在特定代码段中使用新的临时文件夹路径
tempdir(new_tempdir)
# 在此处编写需要使用新临时文件夹的代码


## Pipeline
lavs <- lastz(assemblyTarget, assemblyQuery, 
             outputDir=axtDir,
	      chrsTarget=c(arg[4]),
	      chrsQuery=c(arg[5]),
             distance=c(arg[6]), mc.cores=core)
lavs <- list.files(path=axtDir, pattern="\\.lav$",full.names = TRUE)

psls <- lavToPsl(lavs, removeLav=FALSE, binary="lavToPsl")

chains <- axtChain(psls, assemblyTarget=assemblyTarget,
                   assemblyQuery=assemblyQuery, distance=c(arg[6]),
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

# 恢复默认的临时文件夹路径
tempdir(old_tempdir)
