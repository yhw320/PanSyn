library(syntenet)

setwd("C:/Users/HP/Desktop")
x<-read.table("SynNet-k6s5m25_2cols_infoclusters",header=T,sep="\t")
clusters<-as.data.frame(x)

x1<-read.table("SynNet-k6s5m25_2cols",header=T,sep="\t")
algae_network<-as.data.frame(x1)

cluster_id <- c(7501,7502,7503,11795,11796,11797)

plot_network(algae_network,clusters,cluster_id)
