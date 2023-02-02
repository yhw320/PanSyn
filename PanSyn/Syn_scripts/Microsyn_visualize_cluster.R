arg <- commandArgs(T)
library(tidyverse)
library(magrittr)
library(genoPlotR)
library(tidyverse)
library(magrittr)

setwd(arg[1])
x<-read.table(arg[2],header=T,sep="\t")
data<-as.data.frame(x)
df <- dna_seg(data)

genomes_id <- unique(df$Genome)
genomes_id2 <-sort(genomes_id)
df_dna_seg <- df %>% group_split(Genome) %>% 
  map(~ (.x %>% select(-Genome))) %>% 
  map(~ (.x %>% as.dna_seg())) %>% 
  set_names(nm = genomes_id2)
annot <- lapply(df_dna_seg, function(x) {
  mid <- middle(x)
  annot <- annotation(x1 = mid, text = x$name, rot = 0)
})

ww <- as.numeric(arg[3])
hh <- as.numeric(arg[4])

pdf(file="Conserved_cluster.pdf",width=ww,height=hh)
plot_gene_map(dna_segs = df_dna_seg,dna_seg_label_cex=0.7,annotations = annot,annotation_cex = 0.4,minimum_gap_size = 0.05,scale_cex=0.4)
dev.off()
