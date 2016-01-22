library(flowCore)
library(scales)

ProcessData<-function(address){
  #Stores facs data into a single data frame.
  #Naming convention of fcs files: ID group_factor_level_replicate.fcs (i.e. A01 dh5a_acetoacetate_10_1.fcs).
  setwd(address)
  filenames<-Sys.glob('*.fcs')
  data<-c()
  for (filename in filenames) {
    number=length(exprs(read.FCS(filename))[,3])
    splitname=strsplit(filename,split=' ')[[1]]
    data$id<-c(data$id,rep(splitname[1],each=number))
    data$group<-c(data$group,rep(strsplit(splitname[2],split='_')[[1]][1],each=number))
    data$factor<-c(data$factor,rep(strsplit(splitname[2],split='_')[[1]][2],each=number))
    data$level<-c(data$level,rep(strsplit(splitname[2],split='_')[[1]][3],each=number))
    data$replicate<-c(data$replicate,rep(strsplit(strsplit(splitname[2],split='_')[[1]][4],split='\\.')[[1]][1],each=number))
    data$value<-c(data$value,exprs(read.FCS(filename))[,3])
  }
  data<-data.frame(data$id,data$group,data$factor,data$level,data$replicate,data$value)
  colnames(data)<-c('id','group','factor','level','replicate','value')
  return(data)
}

GetMedians<-function(data){
  #Extracts the median of each fcs file.
  medians<-aggregate(data$value,by=list(data$id,data$replicate,data$level,data$factor,data$group),FUN=median)
  colnames(medians)<-c('id','replicate','level','factor','group','median')
  return(medians)
}

GetMedianMeans<-function(medians,group,factor){
  #Gets the mean and standard deviation of median replicates for a given group and factor.
  subset<-medians[medians$group==group & medians$factor==factor,]
  value<-tapply(subset$median,subset$level,mean)
  dev<-tapply(subset$median,subset$level,sd)
  medianmeans<-as.data.frame(t(rbind(value,dev)))
  medianmeans$level<-as.numeric(rownames(medianmeans))
  return(medianmeans)
}

GetMedianSubset<-function(medians,group,factor){
  #Gets the medians and stdev/sqrt(3) for a group and factor. 
  subset<-medians[medians$group==group & medians$factor==factor,]
  mediansubset<-c()
  value<-subset$median
  level<-as.numeric(as.character(subset$level))
  mediansubset<-data.frame(value,level)
  medianstdev<-tapply(mediansubset$value,mediansubset$level,sd)
  mediansubset$dev<-rep(medianstdev,each=3)/sqrt(3)
  return(mediansubset)
}
