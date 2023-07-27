arg <- commandArgs(T)


library(shape)

x1<-read.table(arg[2],header=F,sep="\t")
hox_data<-as.data.frame(x1)

hox_data[,3] <- as.numeric(hox_data[,3])
hox_data[,2] <- as.numeric(hox_data[,2])

allpos<-rbind(hox_data[,3],hox_data[,2])

xmax <-max(allpos, na.rm = TRUE)
xmin <-min(allpos, na.rm = TRUE)
xmax <- as.numeric(xmax)
xmin <- as.numeric(xmin)

name1<-unique(hox_data[,5])
num<-length(name1)

df <- data.frame(x = c(name1), z = c(1:num))
rownames(df)=df$x

ymin<-1
ymax <-nrow(df)
ymax<-ymax+1

hei <- as.numeric(arg[6])
wid <- as.numeric(arg[7])

# 设置保存文件的目录
directory <- arg[1]

# 设置文件名
filename <- "CNEs.pdf"

# 使用 file.path() 函数构建完整的文件路径
file_path <- file.path(directory, filename)


pdf(file=file_path,height=hei,width=wid)

zheng <- as.numeric(arg[3])
if (zheng==2){
plot(1:5,1:5,xlim=c(-xmax,-xmin),ylim=c(1,ymax),type="n",xlab=arg[4],ylab="",yaxt="n",cex.axis=0.6)
}
if (zheng==1){
plot(1:5,1:5,xlim=c(xmin,xmax),ylim=c(1,ymax),type="n",xlab=arg[4],ylab="",yaxt="n",cex.axis=0.6)
}

#画hox
line_number1<-nrow(hox_data)
while(line_number1>0)
{
spe_name<-hox_data[line_number1,5]
ypos<-df[spe_name,2]
if (zheng==2){
rect(xleft=-c(hox_data[line_number1,2]),ybottom=c(ypos),xright=-c(hox_data[line_number1,3]),ytop=c(ypos+0.6),col=hox_data[line_number1,4],border=hox_data[line_number1,4])
}
if (zheng==1){
rect(xleft=c(hox_data[line_number1,2]),ybottom=c(ypos),xright=c(hox_data[line_number1,3]),ytop=c(ypos+0.6),col=hox_data[line_number1,4],border=hox_data[line_number1,4])
}
line_number1<-line_number1-1
}

#写各个物种的名
rownames(df)<- df[,2]

line_df<-nrow(df)
while(line_df>1)
{
if (zheng==2){
text(x = -xmax,y = line_df+0.8,labels = df[line_df,1],cex=0.9,adj=0)
}
if (zheng==1){
text(x = xmin,y = line_df+0.8,labels = df[line_df,1],cex=0.9,adj=0)
}
line_df<-line_df-1
}
#写hox的名
#rownames(df)<- hox_data[,6]

hox_name1<-read.table(arg[5],header=F,sep="\t")
hox_name<-as.data.frame(hox_name1)


line_df<-nrow(hox_name)
while(line_df>0)
{
xwei=(hox_name[line_df,2] + hox_name[line_df,3])/2
	if (zheng==2){
	text(x = -xwei,y = 1.8,labels = hox_name[line_df,4],cex=0.9)
	}
	if (zheng==1){
	text(x = xwei,y = 1.8,labels = hox_name[line_df,4],cex=0.9)
	}
line_df<-line_df-1
}

dev.off() 
