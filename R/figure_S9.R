source('./process_data.R')
source('./plotting.R')
par(mfrow=c(2,3))
par(mar=c(4.5,5,3,1))

setwd('../data/nissle')
address=getwd()
data<-ProcessData(address)
medians<-GetMedians(data)
groups<-c('nissle','og241','oxb19','patol','loxb19','lpatol')
letters<-c('A)','B)','C)','D)','E)','F)')

count=1
for (item in groups){
  group=item
  if (item=='oxb19' || item=='loxb19'){range<-c(3,6,0,3)}
  else {range<-c(1.5,5,0,3)}
  PlotHist(data,group,'acetoacetate','mM',range,'none',c('black','red'))
  mds=medians[medians$group==group & medians$level==0,]$median
  points(mean(log10(mds)),2.5,col='grey',pch=20,cex=1.5)
  segments(mean(log10(mds))-sd(log10(mds)),2.5,mean(log10(mds))+sd(log10(mds)),2.5,col='grey')
  epsilon<-0.1
  segments(mean(log10(mds))-sd(log10(mds)),2.5-epsilon,mean(log10(mds))-sd(log10(mds)),2.5+epsilon,col='grey')
  segments(mean(log10(mds))+sd(log10(mds)),2.5-epsilon,mean(log10(mds))+sd(log10(mds)),2.5+epsilon,col='grey')
  mds=medians[medians$group==group & medians$level==10,]$median
  points(mean(log10(mds)),2,col='pink',pch=20,cex=1.5)
  segments(mean(log10(mds))-sd(log10(mds)),2,mean(log10(mds))+sd(log10(mds)),2,col='pink')
  epsilon<-0.1
  segments(mean(log10(mds))-sd(log10(mds)),2-epsilon,mean(log10(mds))-sd(log10(mds)),2+epsilon,col='pink')
  segments(mean(log10(mds))+sd(log10(mds)),2-epsilon,mean(log10(mds))+sd(log10(mds)),2+epsilon,col='pink')
  par(new=FALSE)
  par(adj=0)
  title(letters[count],cex.main=1.5)
  par(adj=0.5)
  count=count+1
  if (count==2){legend('topright',1,c('0mM AcAc','10mM AcAc'),col=c('black','red'),box.lty=0,bty='n',lty=c(1,1),lwd=2)}
}

