myFun <- function(data) {
  ListCols <- sapply(data, is.list)
  cbind(data[!ListCols], t(apply(data[ListCols], 1, unlist)))
}#FUNCTION TO UNLIST VARIALBES IN DF


library("lubridate")
library("trajr")
library("sp")

library(dplyr)
library(tibble)

setwd("C:/Users/Rusty/Amazon Drive/MICE")
source("MouseFunctions.R", chdir = F)


IRdir<-"C:/Users/Rusty/Amazon Drive/MICE/Thermal/IRdata"


IRstarttimes<-read.csv("C:/Users/Rusty/Amazon Drive/MICE/Thermal/IRdata/IR_Video_Start_Times.csv")
IRstarttimes<-as.data.frame(IRstarttimes)

IRstarttimes$trial<-sapply(strsplit(as.character(IRstarttimes$file),"_"),"[[",1)#takes 1st element from each stringsplit
IRstarttimes$day<-sapply(strsplit(as.character(IRstarttimes$file),"_"),"[[",3)#takes 3rd element from each stringsplit
IRstarttimes$year<-sapply(strsplit(as.character(IRstarttimes$start_time),"-"),"[[",1)#takes 1st element from each stringsplit
IRstarttimes$time<-sapply(strsplit(as.character(IRstarttimes$start_time)," "),"[[",2)#takes 2nd element from each stringsplit

IRstarttimes$t.hour<-as.numeric(as.character(sapply(strsplit(IRstarttimes$time,":"),"[[",1)))
IRstarttimes$t.min<-as.numeric(as.character(sapply(strsplit(IRstarttimes$time,":"),"[[",2)))
IRstarttimes$t.sec<-as.numeric(as.character(sapply(strsplit(IRstarttimes$time,":"),"[[",3)))



IRstarttimes$dailysecond<-IRstarttimes$t.sec+(IRstarttimes$t.min*60)+(IRstarttimes$t.hour*60*60)



IRfiles<-list.files(IRdir,full.names=TRUE,pattern="shuffle.*\\.csv$")# To date, I've named these megaframes "Summary_xxxxxx.csv"

mice<-c("a1","a2","a3","a4")
IRfileslist<-list()
setwd(IRdir)
for(i in 1:length(IRfiles)){
  coordinfo<-read.csv(IRfiles[[i]])
  colnames(coordinfo)<-NULL
  coordinfo<-myFun(coordinfo)
  coordinfo.labels.a<-as.character(droplevels(unlist(coordinfo[1,])))
  coordinfo.labels.a<-gsub("Mouse","a",coordinfo.labels.a)
  coordinfo.labels.b<-as.character(droplevels(unlist(coordinfo[2,])))
  coordinfo<-coordinfo[-c(1,2),]
  coordinfo<-apply(coordinfo,2,function(x){as.numeric(as.character(x))})
  coordinfo<-as.data.frame(coordinfo)
  
  colnames(coordinfo)<-paste(as.character(coordinfo.labels.a),as.character(coordinfo.labels.b),sep='.')
  colnames(coordinfo)[1]<-"frame"
  ###########################################
  #time-match info
  filename<-strsplit(IRfiles[[i]],"/")[[1]][8]
  trial<-strsplit(filename,"_")[[1]][1]
  day<-gsub("DLC","",strsplit(filename,"_")[[1]][3])
  timeinfo<-IRstarttimes[which(IRstarttimes$trial==trial & IRstarttimes$day==day),]
  ##########################
  for(m in 1:4){
    tempinfo<-coordinfo[,grep(mice[[m]],colnames(coordinfo))]
    tempinfo[,1]<-ifelse(tempinfo[,3]>0.9,tempinfo[,1],NA)
    tempinfo[,2]<-ifelse(tempinfo[,3]>0.9,tempinfo[,2],NA)
    #IMPUTE MISSING X/Y COORDINATES FOR MOUSE LOCATIONS
    tempinfo[,c(1,2)]<-apply(tempinfo[,c(1, 2)],2,
                             function(X){approxfun(seq_along(X),X)(seq_along(X))})
    
    tempinfo<-cbind(coordinfo$frame,tempinfo[,c(1,2)])
    colnames(tempinfo)[1]<-"frame"
    # trj <- TrajFromCoords(track=tempinfo[,c(2,3,1)],
    #                         #xCol="a1x0",yCol="a1y0",
    #                         timeCol="frame",spatial="pixels")
    # 
    
    if(m==1){
      micey<-tempinfo
    } else {
      micey<-cbind(micey,tempinfo[,c(2,3)]) 
    }
  }
  
  #calculats pairwise dist between each pair of mice at every instant
  CombinedInfo<-pairwise.Distances(micey,inds=4,pix.cm=4.2924)#Custom function located in MouseFunctions.R
  
  #are a given pair of mice associating at a given time (distance < AssociatingDistance)
  CombinedInfo<-Associate.Identifier(CombinedInfo,AssociatingDistance=10)#Custom function located in MouseFunctions.R
  
  #
  #Creates large, mouse-to-mouse dataset at each step ()
  pdf(paste(trial,day,"MouseAngles.pdf",sep="_"),width= 8, height= 8,family="NimbusRom")
  par(mfrow=c(4,4))# rows, columns
    MouseAngles<-Mouse2Mouse(xyvalues.allpairs=CombinedInfo[,c(2:9)],
                             pairwisedistances=CombinedInfo[,grep("dist",colnames(CombinedInfo))],
                             n.inds=4,shrinksize=.85)
  dev.off()
  
  #############
  
  Big<-cbind(CombinedInfo[-1,],MouseAngles)
  #Categorizes behavior
  behaviors<-BehaviorCategorizer(Big,approach.angle=40,
                                 leave.angle=90,integratesteps=10,n.inds=4,
                                 walk.speed=10,stationary.noise=.5)
  HugeMouse<-cbind(Big,behaviors)
  #colnames(HugeMouse)
  #summary(HugeMouse)
  
  es<-(timeinfo$dailysecond):(36000+timeinfo$dailysecond) #make it up to 10 hrs longer (36000 sec)
  ds<-rep(es,each=10) # 10x of each dailysecond
  ds<-ds[c(1:nrow(HugeMouse))] # make it only as long as the HugeMouse
  HugeMouse<-cbind(ds,HugeMouse)
  
  IRbehavior.for.matching.thermal<-HugeMouseDownSampler(HugeMouse)
 
  IRfileslist[[i]]<-IRbehavior.for.matching.thermal
  names(IRfileslist)[i]<-paste(trial,day,"IR",sep="_") 
  
  save(IRfileslist,file="IRfileslist.Rdata")
}


save(IRfileslist,file="IRfileslist.Rdata")








