









library(sm)
library(XLConnect)

#?readWorksheetFromFile
setwd("C:/Users/Rusty/Amazon Drive/MICE")

mousedata<-readWorksheetFromFile("8x8_Master.xlsx",sheet="Single Quad Experiments",useCachedValues=TRUE)


mousedata$id<-factor(mousedata$id)
mousedata$X8x8_age<-as.numeric(as.character(mousedata$X8x8_age))
mousedata$X8x8_pre_mass.g.<-as.numeric(as.character(mousedata$X8x8_pre_mass.g.))
mousedata$X8x8_post_mass.g.<-as.numeric(as.character(mousedata$X8x8_post_mass.g.))
mousedata$AGD.mm.<-as.numeric(as.character(mousedata$AGD.mm.))
mousedata$gonad_mass.g.<-as.numeric(as.character(mousedata$gonad_mass.g.))
mousedata$OFT_age<-as.numeric(as.character(mousedata$OFT_age))
mousedata$OFT_mass.g.<-as.numeric(as.character(mousedata$OFT_mass.g.))






strains<-unique(mousedata$strain)
behaviorstoplot<-colnames(mousedata)[c(31,37,33,29)]


pdf("OFT_mousebehavior.pdf",width= 24, height= 24,family="NimbusRom")
  par(mfrow=c(5,4))# rows, columns

for(cv in 1:length(strains)){
  
  dis.strain<-strains[cv]
  
  thisstrain<-subset(mousedata,strain==dis.strain)
  
  for(cc in 1:length(behaviorstoplot)){
    
    thisstrain.behavior<-thisstrain[,c("sex",behaviorstoplot[cc])]
    thisstrain.behavior<-na.omit(thisstrain.behavior)
    
    
    if(nrow(thisstrain.behavior)>2){
    
    # create value labels 
    sex.f <- factor(thisstrain.behavior$sex, levels= c("f","m"),
                    labels = c("female", "male")) 
    
    # plot densities 
    sm.density.compare(thisstrain.behavior[,2], sex.f, xlab=behaviorstoplot[cc],col=c("blue","purple"),lwd=3)
    title(main=paste(dis.strain,behaviorstoplot[cc],sep=":"))
    
    # add legend via mouse click
    colfill<-c(2:(2+length(levels(sex.f)))) 
    legend("topleft", levels(sex.f), col=c("blue","purple"),lwd=3,lty=c(1,5),seg.len=4)
    
    }
    
  }
  
}

dev.off()





