setwd('../data')
RFU<-read.csv('./dh5a_nissle.csv',sep=',',header=TRUE)[2:101]
OD<-read.csv('./dh5a_nissle.csv',sep=',',header=TRUE)[102:201]
NormRFU<-RFU/OD
times<-t(read.csv('./dh5a_nissle.csv',sep=',',header=FALSE)[1,2:101])
names<-as.character(read.csv('./dh5a_nissle.csv',sep=',',header=TRUE)[,1])
NormRFU<-data.frame(t(NormRFU))
colnames(NormRFU)<-c(names)

ylim1=0.9
ylim2=1.2
par(mar=c(4.5,5,4,1))
par(mfrow=c(2,2))
for (index in c(1,7,13,19)){
  y1<-apply(NormRFU[,index:(index+2)],FUN=mean,MARGIN=1)
  y2<-apply(NormRFU[,(index+3):(index+5)],FUN=mean,MARGIN=1)
  y<-y2/y1
  sd1<-apply(NormRFU[,index:(index+2)],FUN=sd,MARGIN=1)
  sd2<-apply(NormRFU[,(index+3):(index+5)],FUN=sd,MARGIN=1)
  sd<-sqrt((sd1/y1)^2+(sd2/y2)^2)*(y2/y1)
  x<-times
  plot(x,y,pch=1,col='white',xlim=c(0,990),ylim=c(ylim1,ylim2),ylab='Norm. RFU Ratio',xlab='time after induction (mins)',cex=0.8,cex.lab=1.5,cex.axis=1.3)
  grid(NULL,NULL,lwd=1)
  i=1
  for (item in sd){
    segments(x[i],y[i]-item,x[i],y[i]+item,col='black')
    segments(x[i]-5,y[i]-item,x[i]+5,y[i]-item,col='black')
    segments(x[i]-5,y[i]+item,x[i]+5,y[i]+item,col='black')
    i=i+1
  }
  par(new=TRUE)
  plot(x,y,pch=21,col='black',bg='white',xlim=c(0,990),ylim=c(ylim1,ylim2),xaxt='n',yaxt='n',bty='n',xlab='',ylab='',main='',cex=0.8)
  par(adj=0)
  if (index==1){
    title('A)',cex.main=1.7)
    ylim1=0.9
    ylim2=1.2
  }
  if (index==7){
    title('B)',cex.main=1.7)
    ylim1=0.9
    ylim2=1.5
  }
  if (index==13){
    title('C)',cex.main=1.7)
  }
  if (index==19){
    title('D)',cex.main=1.7)
  }
  par(adj=0.5)
}
