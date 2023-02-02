arg <- commandArgs(T)
setwd(arg[1])
sink("Pvalue-result.txt")#输出的文件名
x<-read.table(arg[2],header=T,sep="\t")
data<-as.data.frame(x)

#统计列
num_col<-ncol(data)
#统计行
num_row<-nrow(data)

i_row=1

while(i_row<num_row)
{
	i_col=2
	while(i_col<num_col){
		a<-data[i_row,i_col]
		b<-data[i_row,num_col]-a
		c<-data[num_row,i_col]-a
		d<-data[num_row,num_col]-a-b-c
		y1<-c(a,b,c,d)
		alle<-matrix(y1, nrow=2)
		a<-fisher.test(alle,alternative ="greater")
		if (i_col==num_col-1){
			cat(a$p.value,'\n')
		}
		else{
			cat(a$p.value,'\t')
		}
		i_col=i_col+1
	}
	i_row=i_row+1
}
sink()
