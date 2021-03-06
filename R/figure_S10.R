source('./process_data.R')
source('./plotting.R')
par(mfrow=c(1,4))
par(mar=c(4.5,5,3,1))

setwd('../data/scfa')
address=getwd()
data<-ProcessData(address)
medians<-GetMedians(data)

range<-c(2,4.5,0,2.5)

PlotHist(data,'patol','acetoacetate','mM',range,'none',c('grey','red'))
legend('topright',1,c('0mM','10mM'),col=c('grey','red'),box.lty=0,bty='n',lty=c(1,1),lwd=2)
par(adj=0)
title('A)',cex.main=1.5)
par(adj=0.5)
par(new=FALSE)

PlotHist(data,'patol','acetate','mM',range,'none',c('grey','orange','red'))
legend('topright',1,c('0mM','10mM','50mM'),col=c('grey','orange','red'),box.lty=0,bty='n',lty=c(1,1),lwd=2)
par(adj=0)
title('B)',cex.main=1.5)
par(adj=0.5)
par(new=FALSE)

PlotHist(data,'patol','proprionate','mM',range,'none',c('grey','orange','red'))
legend('topright',1,c('0mM','10mM','50mM'),col=c('grey','orange','red'),box.lty=0,bty='n',lty=c(1,1),lwd=2)
par(adj=0)
title('C)',cex.main=1.5)
par(adj=0.5)
par(new=FALSE)

PlotHist(data,'patol','butyrate','mM',range,'none',c('grey','orange','red'))
legend('topright',1,c('0mM','10mM','50mM'),col=c('grey','orange','red'),box.lty=0,bty='n',lty=c(1,1),lwd=2)
par(adj=0)
title('D)',cex.main=1.5)
par(adj=0.5)
par(new=FALSE)
