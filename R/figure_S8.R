source('./process_data.R')
source('./plotting.R')
par(mfrow=c(3,4))
par(mar=c(4.5,5,3,1))

setwd('../data/facs')
address=getwd()
data<-ProcessData(address)
medians<-GetMedians(data)
groups<-c('dh5a','og241','oxb19','patog','patok','patol','patom')
letters<-c('A)','B)','C)','D)','E)','F)','G)')

count=1
for (item in groups){
  group=item
  if (item=='oxb19'){range<-c(4,7.5,0,2.5)}
  else {range<-c(2,5,0,2.5)}
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
  if (count==2){legend('topright',1,c('0mM','10mM'),col=c('black','red'),box.lty=0,bty='n',lty=c(1,1),lwd=2)}
}

setwd('../facs_low')
address=getwd()
data<-ProcessData(address)
medians<-GetMedians(data)
groups<-c('log241','loxb19','lpatok','lpatol','lpatom')
letters<-c('H)','I)','J)','K)','L)')

count=1
for (item in groups){
  group=item
  range<-c(2,5,0,2.5)
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
}

