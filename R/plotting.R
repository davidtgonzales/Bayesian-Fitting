PlotHist<-function(data,group,factor,units,range,side,colors){
  #Plots histograms of each level in a group and factor in the same plot.
  subset<-data[data$group==group & data$factor==factor,]
  #colors<-rev(heat.colors(length(unique(subset$level))))
  xrange<-c(range[1],range[2])
  yrange<-c(range[3],range[4])
  border<-c(1)
  color<-c(1)
  levels<-sort(unique(levels(subset$level)))
  for (i in levels){
    for (j in unique(subset$replicate)[1:length(unique(subset$replicate))]){
      density<-density(log10(subset[subset$level==i & subset$replicate==j,]$value),log='x',na.rm=TRUE)
      if (border==1){
        plot(density,main='',xlab='log10(RFU)',ylab='Density',lwd=2,ylim=c(yrange[1],yrange[2]),xlim=c(xrange[1],xrange[2]),col=colors[color],cex.lab=1.7,cex.axis=1.5,cex=1.3)
        grid(NULL,NULL,lwd=1)
        border<-c(2)
        par(new=TRUE)
      }
      if (border!=1){
        plot(density,xaxt='n',yaxt='n',bty='n',xlab='',ylab='',main='',lwd=2,ylim=c(yrange[1],yrange[2]),xlim=c(xrange[1],xrange[2]),col=colors[color])
        par(new=TRUE)
      }
    }
    par(new=TRUE)
    color<-c(color+1)
  }
  if (side=='right'){
    legend('topright',yrange[2],paste(levels,factor,sep=paste(units,' ',sep='')),col=colors,box.lty=0,bty='n',lty=1,lwd=3)
  }
  else if (side=='left'){
    legend('topleft',yrange[2],paste(levels,factor,sep=paste(units,' ',sep='')),col=colors,box.lty=0,bty='n',lty=1,lwd=3)
  }
}

library(scales)
PlotCurve<-function(data,times,color,limits,borders,xlabel,ylabel,mainlabel){
  lightcolor<-alpha(color,1)
  if (borders=='no'){
    par(new=TRUE)
    plot(rep(times,3),t(data),pch=1,xlim=c(limits[1],limits[2]),ylim=c(limits[3],limits[4]),col=lightcolor,xaxt='n',yaxt='n',bty='n',xlab='',ylab='',main='',cex=0.8)
    par(new=FALSE)
  }
  else if (borders=='yes'){
    plot(rep(times,3),t(data),pch=1,xlim=c(limits[1],limits[2]),ylim=c(limits[3],limits[4]),col=lightcolor,xlab=xlabel,ylab=ylabel,main='',cex.lab=1.7,cex.axis=1.5,cex=0.8)
  }
  grid(NULL, NULL, lwd = 1)
}
