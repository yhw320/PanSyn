arg <- commandArgs(T)


output_dir <- arg[1]  # 获取命令行参数中的输出文件目录
file_name <- "Pvalue-result.txt"  # 输出文件的文件名

output_file <- file.path(output_dir, file_name)  # 生成完整的输出文件路径

sink(output_file)  # 设置输出文件路径
x <- read.table(arg[2], header = TRUE, sep = "\t")
data <- as.data.frame(x)



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
