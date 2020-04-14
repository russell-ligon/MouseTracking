


library(gdata)
library(zoo)
library(circular)


####DOWNSAMPLES DATA to 1HTZ
onesecond.downsampler<-function(datatodownsample=profileplainpee,startingnames=profileplainpee$label){
  
  #print(paste("Downsampling ",dir,sep=''))
  #frameinfo<-lapply(FullFiles.FullNames,function(x) strsplit(x,"2018-"))

  if(length(grep("2019",startingnames[1]))==0){
    frameinfo<-strsplit(startingnames,"2019-")
  } else {
    frameinfo<-strsplit(startingnames,"2018-")
  }
  
  
  alltimes<-lapply(frameinfo,function(x) strsplit(x[2],".csv"))
  oknow<-lapply(alltimes,function(x) strsplit(x[[1]],"_"))
  oknow2<-unlist(lapply(oknow,function(x) x[[1]][c(2)]))
  oknow3<-unlist(lapply(oknow,function(x) x[[1]][c(1)]))
  
  fileinformation<-data.frame(matrix(unlist(alltimes),nrow=length(alltimes),byrow=T));colnames(fileinformation)<-"V1"
  fileinformation$day<-oknow3
  fileinformation$tosecond<-oknow2
  fileinformation$check<-paste(oknow3,oknow2,sep="-")
  fileinformation$unique.time<-!duplicated(fileinformation$check)      
  
  
  #deleters<-fileinformation[which(!fileinformation$unique.time),]
  
  downsampled<-datatodownsample[which(fileinformation$unique.time),]
  return(downsampled)
}
#downsample just list of names
onesecondgroups<-function(startingnames){
  if(length(grep("2019",startingnames[1]))==0){
    frameinfo<-strsplit(startingnames,"2018-")
  } else {
    frameinfo<-strsplit(startingnames,"2019-")
  }
  
  alltimes<-lapply(frameinfo,function(x) strsplit(x[2],".csv"))
  oknow<-lapply(alltimes,function(x) strsplit(x[[1]],"_"))
  oknow2<-unlist(lapply(oknow,function(x) x[[1]][c(2)]))
  oknow3<-unlist(lapply(oknow,function(x) x[[1]][c(1)]))
  
  fileinformation<-data.frame(matrix(unlist(alltimes),nrow=length(alltimes),byrow=T));colnames(fileinformation)<-"V1"
  fileinformation$day<-oknow3
  fileinformation$tosecond<-oknow2
  fileinformation$check<-paste(oknow3,oknow2,sep="-")
  return(fileinformation$check)
}
#just return downsample indices
onesecond.index<-function(startingnames){
  if(length(grep("2019",startingnames[1]))==0){
    frameinfo<-strsplit(startingnames,"2018-")
  } else {
    frameinfo<-strsplit(startingnames,"2019-")
  }
  
  alltimes<-lapply(frameinfo,function(x) strsplit(x[2],".csv"))
  oknow<-lapply(alltimes,function(x) strsplit(x[[1]],"_"))
  oknow2<-unlist(lapply(oknow,function(x) x[[1]][c(2)]))
  oknow3<-unlist(lapply(oknow,function(x) x[[1]][c(1)]))
  
  fileinformation<-data.frame(matrix(unlist(alltimes),nrow=length(alltimes),byrow=T));colnames(fileinformation)<-"V1"
  fileinformation$day<-oknow3
  fileinformation$tosecond<-oknow2
  fileinformation$check<-paste(oknow3,oknow2,sep="-")
  fileinformation$unique.time<-!duplicated(fileinformation$check) 
  return(which(fileinformation$unique.time))
  
}

#takes dataframe (use keeplongA13, keeplongA24), and trial/day specific roi information to add the roi where pee happens
#RETURNS [[1]] dataframe with ROI-identified pee coords & [[2]] dataframe with info on errors
peeroi<-function(dfwp,inforoi){
  library(sp)
  dfwp$peeroi<-NA
  dfwp$peeroiFULL<-NA
  dfhold<-dfwp
  
  dfwp2<-dfwp[which(!is.na(dfwp$mTemp)),]
  replacementrowvalue<-which(!is.na(dfwp$mTemp))
  
  #pulls just the roi info for the correct camera
  if(length(grep("a1",colnames(dfwp)))>1){
    inforoi2<-inforoi[which(inforoi$arena=="a1" | inforoi$arena=="a3"),]
  }
  if(length(grep("a2",colnames(dfwp)))>1){
    inforoi2<-inforoi[which(inforoi$arena=="a2" | inforoi$arena=="a4"),]
  }
  
  
  if(nrow(dfwp2)==0){
    dfhold<-dfhold
  } else {
    errorflag<-1
    for(r in 1:nrow(dfwp2)){
      ne<-0
      for(s in 1:nrow(inforoi2)){
        
        e<-point.in.polygon(dfwp2$true.x[r],dfwp2$true.y[r],pol.x=c(inforoi2[s,c(6,8,10,12)]),pol.y=c(inforoi2[s,c(7,9,11,13)]))
        #e2<-ifelse(e==1 & dfwp2[r,"quad"]==inforoi[s,"arena"],1,0)
        ne<-c(ne,e)
      }
      ne<-ne[-1] #drop first entry, which was created as 0 before s loop
      roispossible<-inforoi2$ROI_name[ne==1]
      roispossible<-unique(roispossible)
      fullpossiblerois<-paste(as.character(roispossible),collapse=".")
      if(length(roispossible)>1){#possible for a) waters, which are also on walls, and water is given precedence
        # and b) corners, and corners given precedence
        #print(as.character(roispossible))
        if(length(roispossible[grep("water",roispossible)])>0){
          roispossible<-roispossible[grep("water",roispossible)]
        }
        if(length(roispossible[grep("barrier",roispossible)])==2){
          roinam<-as.character(roispossible[1])
          roispossible<-paste0(strsplit(roinam,"_")[[1]][1],"_central_corner")
        }
        if(length(roispossible[grep("corner",roispossible)])>0){
          roispossible<-roispossible[grep("corner",roispossible)]
        }
        if(((length(roispossible[grep("barrier",roispossible)])>0)+(length(roispossible[grep("wall",roispossible)])))==2){
          roispossible<-roispossible[grep("barrier",roispossible)]
        }
        
      }
      
      #if multiple 'corner' rois possible, use the one that matches the quad identified by me, earlier in the code
      if(length(roispossible[grep("corner",roispossible)])>1){
        roispossible<-roispossible[grep(dfwp2[r,"quad"],roispossible)]
        #if the above trim doesn't work, length will still be >1
        #IN which case we take the lead quad
        if(length(roispossible)>1){
          dropindex<-which(roispossible==roispossible[grep(paste0("_",dfwp2[r,"quad"]),roispossible)])
          roispossible<-roispossible[-dropindex]
        }
      }
      
      if(length(roispossible)==0){
        roispossible<-"generalmid"
      }
      if(length(fullpossiblerois)==0){
        fullpossiblerois<-"NOMATCH"
      }
      
      keepnewrow<-dfwp2[r,]
      keepnewrow$peeroi<-as.character(roispossible)
      keepnewrow$peeroiFULL<-fullpossiblerois
      replacindex<-replacementrowvalue[r]
      
      dfhold[replacindex,]<-keepnewrow
      
      if(r==1){
        keepnewframe<-keepnewrow
      } else {
        keepnewframe<-rbind(keepnewframe,keepnewrow)
      }
      
      #capture roi and quad mismatch instances, but ignore instances where the inferred location is generalmid (impossible to compare)
      if(roispossible!="generalmid"){
        if((dfwp2[r,"quad"]!=strsplit(as.character(roispossible),"_")[[1]][1])){
          errorstring<-data.frame(c(inforoi2[1,c("trial","camera","day")],as.character(dfwp2[r,"quad"]),as.character(roispossible),as.character(fullpossiblerois)))
          colnames(errorstring)[c(4:6)]<-c("ralQuad","roiID","fullmatches")
          #print(r)
          if(errorflag==1){
            allerrors<-errorstring
            errorflag<-errorflag+1
          } else {
            allerrors<-rbind(allerrors,errorstring)
          }
        }#close if testing whether the first quad component of the roi name is the same as the RAL identified pee quad
        
      }#close if testing whether 'generalmid'
      
    }
    
    dfhold$peeroi<-factor(dfhold$peeroi)
    #summary(dfhold)
  }
  
  if(!exists("allerrors")){
    errorstring<-data.frame(c(inforoi[1,c("trial","camera","day")],as.character(dfwp2[1,"quad"]),as.character("moop")))
    colnames(errorstring)[c(4:5)]<-c("ralQuad","roiID")
    allerrors<-errorstring[0,]
  }
  
  wantedlist<-list(dfhold,allerrors)
  
  return(wantedlist)
  
  
}

######renametherois
# function that takes old rois and converts them, along with relevant infor from 'thistrial' to create new universal ROI labels
renametherois<-function(v.of.rois,thistrial,sex.of.focal){
  # v.of.rois  = vector of rois in original format, eg. a4_a2_barrier
  # thistrial = subset of FullMegaInfo for this trial
  # sex.of.focal = simply the sex (f or m) of the focal mouse in a given quad
  
  v2<-gsub("-","_",v.of.rois)
  
  hooke<-strsplit(as.character(v2),"_")
  neighbor.vec<-unlist(lapply(hooke,function(x){
    if(length(x)==3){
      neighbor<-x[2]
    } else {
      neighbor<-""
    }}))
  
  roitype.vec<-unlist(lapply(hooke,function(x){x[length(x)] }))
  
  neighbor.sex<-unlist(lapply(neighbor.vec,function(x){
    if(x!=""){
      if(any(thistrial$quad==x)){
        bing<-as.character(thistrial[which(thistrial$quad==x),"sex"])
      } else {
        bing<-""
      }
      
    } else {
      bing<-""
    }
    return(bing)
  }
  )
  )
  
  sexcompare<-data.frame(cbind(as.character(sex.of.focal),as.character(neighbor.sex)))
  colnames(sexcompare)<-c("focalsex","compsex")
  sexcompare$focalsex<-as.character(sexcompare$focalsex)
  sexcompare$compsex<-as.character(sexcompare$compsex)
  
  sexcompare$sexcat<-ifelse(sexcompare$compsex=='','',
                            ifelse(sexcompare$focalsex==sexcompare$compsex,"SS","OS"))
  
  neightbor2<-neighbor.vec
  neightbor2[!is.na(match(neighbor.vec,c("a1","a2","a3","a4")))]<-""
  addmodifier<-neightbor2
  
  peeishappeninghere<-paste0(sexcompare$sexcat,addmodifier,roitype.vec)
  
  return(peeishappeninghere)
  
}


#Barplotting function, with error bars
barplotwitherrorbars<-function(dataframe,valuecolumn,groupingvariablecolumn,secondarygroupingvariable=NA,plottingcolors,cexName=1.3){
  library(dplyr)
  #dataframe = dataframe to summarize, multiple observations = rows, variables = columns
  
  #valuecolumn = name of column from a dataframe that you want plotted, e.g. "height"
  
  #groupingvariablecolumn (e.g. "sex")
  
  #plottingcolors = vector of color values corresponding to alphabetic ording of grouping variables
  
  #secondarygroupingvariable = if you want a double split (e.g. Sex X Strains), this is where you define that variable
  
  if(!exists("cexName")){cexName=2}
  
  if(!is.na(secondarygroupingvariable)){
    possiblecombos<-unique(expand.grid(dataframe[,c(groupingvariablecolumn,secondarygroupingvariable)]))
    
    flag<-1
    
    ns<-list()
    sems<-list()
    means<-list()
    varnams<-list()
    
    types.primary<-unique(dataframe[,groupingvariablecolumn])
    types.secondary<-unique(dataframe[,secondarygroupingvariable])
    
    for(r in 1:length(types.secondary)){
      secondaryval<-types.secondary[r]
      subseconddata<-dataframe[which(dataframe[,secondarygroupingvariable]==secondaryval),]
      
      for(s in 1:length(types.primary)){
        primaryval<-types.primary[s]
        subprim<-subseconddata[which(subseconddata[,groupingvariablecolumn]==primaryval),]
        
        
        
        
        
        ns[[flag]]<-nrow(subprim)
        varnams[[flag]]<-paste(possiblecombos[flag,1],possiblecombos[flag,2],sep=".")
        calcvalue<-na.omit(subprim[,valuecolumn])
        #computation of the standard error of the mean
        sems[[flag]]<-sd(calcvalue)/sqrt(length(calcvalue))
        means[[flag]]<-mean(calcvalue)
        
        
        
        flag<-flag+1
        
      }
      
    }
    
    standardErrors <- unlist(sems)
    means<-unlist(means)
    names<-unlist(varnams)
    samplesizes<-unlist(ns)
    
    plotTop <- max(means+standardErrors*2)
    barCenters <- barplot(means, names.arg=names, col=plottingcolors, las=1, ylim=c(0,plotTop),ylab=valuecolumn,cex.names = cexName,las=2)
    segments(barCenters, means-standardErrors*2, barCenters, means+standardErrors*2, lwd=2)
    text(x = barCenters+.2, y= means-means*.25, label = samplesizes, pos = 3, cex = 0.8, col = "blue")
    
  } else {
    #default, bar plot with single grouping variable and SEM bars
    types<-dataframe[,groupingvariablecolumn]
    
    ns<-list()
    sems<-list()
    means<-list()
    varnams<-list()
    for(g in 1:nlevels(types)){
      g2<-levels(types)[g]
      varnams[[g]]<-g2
      calcvalue<-na.omit(dataframe[which(dataframe[,groupingvariablecolumn]==g2),valuecolumn])
      #computation of the standard error of the mean
      sems[[g]]<-sd(calcvalue)/sqrt(length(calcvalue))
      means[[g]]<-mean(calcvalue)
      ns[[g]]<-length(calcvalue)
    }
    
    standardErrors <- unlist(sems)
    means<-unlist(means)
    names<-unlist(varnams)
    samplesizes<-unlist(ns)
    
    plotTop <- max(means+standardErrors*2)
    barCenters <- barplot(means, names.arg=names, col=plottingcolors, las=1, ylim=c(0,plotTop),ylab=valuecolumn,cex.names = cexName)
    segments(barCenters, means-standardErrors*2, barCenters, means+standardErrors*2, lwd=2)
    text(x = barCenters+.2, y= means-means*.25, label = samplesizes, pos = 3, cex = 0.8, col = "blue")
  }
}



# CreateCompositeList -----------------------------------------------------
#Loops through each folder in 'directories'
# pulls files corresponding to tracking informaiton (location.names)
# and ROI information (Roifiles)
# then names the information from these files and amalgamates into a data.frame called CombinedInfo
# then puts that dataframe into the list AllFoldersList at position zz

CreateCompositeList<-function(directories,pix.cm=4.2924,AssociatingDistance.cm=10,
                                xytrackpattern="fixed.csv",roipattern="-region",roundxy=TRUE){
  #directories should include all folders/subfolder with tracking data
  #pix.cm defines the pixels/cm correction
  #AssociatingDistance.cm Defines cm distance that corresponding to an 'association'
  AllFoldersList<-list() #Creates empty list which will get filled iteratively with loop
  flag<-1 #sets flag, which will only increase when looping through folders containing the "correct" kind of data
  
  
  for (zz in 1:length(directories)){
    FolderInfo1<-directories[zz]#pull full 
    FolderInfo<-strsplit(FolderInfo1,"/")[[1]][length(strsplit(FolderInfo1,"/")[[1]])]
    FolderInfo<-paste(strsplit(FolderInfo1,"/")[[1]][(length(strsplit(FolderInfo1,"/")[[1]])-1)],FolderInfo,sep='.')
    location.names<-list.files(directories[zz],full.names=TRUE,pattern=xytrackpattern)#change pattern to only pull location csvs
    
    if(length(location.names)>0){ #Only if in the right kind of folder, i.e. containing ...fixed.csv files, run the rest, otherwise, skip
      
      Roifiles<-list.files(directories[zz],full.names = TRUE,pattern=roipattern)#refine to pull roi csvs
      Roifilenames<-list.files(directories[zz],full.names = FALSE,pattern=roipattern)#refine to pull roi csvs
      individuals<-length(location.names)
      
      for(Caleb in 1:individuals){
        Location<-read.csv(location.names[Caleb])
        colnames(Location)[2:3]<-paste(colnames(Location)[2:3],".A",Caleb,sep='')
        Location[,3]<-(Location[,3]-1080)*(-1)
        
        if(Caleb>1){
          CombinedInfo<-merge(CombinedInfo,Location,by="position")
        } else {
          CombinedInfo<-Location
        }
        
      }
      
      
      #CombinedInfo[,c(2:ncol(CombinedInfo))]<-CombinedInfo[,c(2:ncol(CombinedInfo))]/(pix.cm)
      if(roundxy==TRUE){
      CombinedInfo<-round(CombinedInfo)
      }
        if(individuals>1){
          CombinedInfo<-pairwise.Distances(CombinedInfo,individuals)#Custom function located in MouseFunctions.R
          ncomparisons<-2*(individuals-1) #Calculates number of unique dyadic comparisons based on the n of individuals
          dister<-(individuals*2+2)
          CombinedInfo[,c(dister:(dister+ncomparisons-1))]<-(CombinedInfo[,c(dister:(dister+ncomparisons-1))])/(pix.cm)
          CombinedInfo<-Associate.Identifier(CombinedInfo,AssociatingDistance.cm)#Custom function located in MouseFunctions.R
        }
      if(length(Roifilenames)>1){
        for(Caitlin in 1:length(Roifiles)){
          ROI<-read.csv(Roifiles[Caitlin])
          roiname<-Roifilenames[Caitlin]
          colnames(ROI)[2]<-roiname
          CombinedInfo<-merge(CombinedInfo,ROI, by="position")
        }
      }
      
      AllFoldersList[[flag]]<-CombinedInfo #puts CombinedInfo dataframe into AllFoldersList at position zz
      names(AllFoldersList)[[flag]]<-FolderInfo #applies name of folder to list element
      flag<-flag+1
    }
  }  
  
  return(AllFoldersList)
}



# BehaviorCategorizer -----------------------------------------------------
#Big =
# approach.angle = acceptable width of angle for approaching another mouse
# leave.angle = acceptable width of angle for leaving another mouse
# integratesteps = n of frames over which to calculate rolling averages for angles and delta pairwise distances for calculating approaching/leaving
# walk.speed = cm/sec threshold, above which = running

BehaviorCategorizer<-function(Big,approach.angle=90,leave.angle=90,integratesteps=10,n.inds=4,walk.speed=10,stationary.noise=.5){
  orderedcomparisons<-2*(n.inds-1)*2
  OrderedMicePairs<-colnames(Big)[grep(":",colnames(Big))]
  BigSub<-Big[,c(grep(":",colnames(Big)),grep("Delta",colnames(Big)))]
  Bog<-rollapply(BigSub,width=integratesteps,mean,na.rm=TRUE,partial=TRUE,align="left")
  MouseDisp<-Big[,grep("step",colnames(Big))]
  RunMouse<-rollapply(MouseDisp,width=integratesteps,sum,na.rm=TRUE,partial=TRUE,align="left")
  colnames(RunMouse)<-gsub("step","cm/s",colnames(RunMouse))  
  #left alignment means that the locations we identify for behaviors (e.g. approaches)
  # will corespond to the 'start' of the behavior
  
  
  movingbehavior<-apply(RunMouse,2,function(x) ifelse(x<stationary.noise,"stationary",
                                                      ifelse(x>walk.speed,"running","walking")))
  colnames(movingbehavior)<-gsub("cm/s","movement",colnames(movingbehavior))
  
  #Approaching classifier
  top.angle<-approach.angle/2
  bottom.angle<-(approach.angle/2)*-1
  ApproachAngleClassifier<-function(angles){
    angleapproach<-ifelse(angles>bottom.angle & angles<top.angle,1,0)
  }
  
  #Leaving classifier
  Ltop.angle<-180-(leave.angle/2)
  Lbottom.angle<-(-180)+(leave.angle/2)
  LeaveAngleClassifier<-function(angles){
    anglealeave<-ifelse(angles>Ltop.angle | angles<Lbottom.angle,1,0)
  }
  
  Approaching<-apply(Bog[,c(grep(":",colnames(Bog)))],2,ApproachAngleClassifier)
  colnames(Approaching)<-gsub("-angle","pointedtoward",colnames(Approaching))
  Leaving<-apply(Bog[,c(grep(":",colnames(Bog)))],2,LeaveAngleClassifier)
  colnames(Leaving)<-gsub("-angle","pointedaway",colnames(Leaving))
  
  socialmovements<-cbind(Approaching, Leaving,Bog[,c(grep(":",colnames(Bog)))],Bog[,c(grep("Delta",colnames(Bog)))])
  
  
  oc<-strsplit(OrderedMicePairs,"-")
  oc2<-lapply(oc,function(x) strsplit(x[1],":"))
  oc3<-lapply(oc2,function(x) unlist(x))
  
  
  for(rory in 1:orderedcomparisons){
    plusr<-rory+orderedcomparisons
    
    approachingthismouse<-ifelse(socialmovements[,rory]==1 & 
                                   socialmovements[,c(grepl(oc3[[rory]][1],colnames(socialmovements))&
                                                                                  grepl(oc3[[rory]][2],colnames(socialmovements))&
                                                                                  grepl("Delta",colnames(socialmovements)))]<0 & 
                                                                                  movingbehavior[,grepl(oc3[[rory]][1],colnames(movingbehavior))]!="stationary",1,0)
    approachingthismouse<-as.data.frame(approachingthismouse,ncol=1)
    colnames(approachingthismouse)<-  gsub(":",".app.",gsub("-angle","",OrderedMicePairs))[rory]

    
    if(rory<2){
      approachingbehavior<-approachingthismouse
    } else {
      approachingbehavior<-cbind(approachingbehavior,approachingthismouse)
    }
 
    leavingthismouse<-ifelse(socialmovements[,plusr]==1 & socialmovements[,c(grepl(oc3[[rory]][1],colnames(socialmovements))&
                                                                               grepl(oc3[[rory]][2],colnames(socialmovements))&
                                                                               grepl("Delta",colnames(socialmovements)))]>0 &
                                                                               movingbehavior[,grepl(oc3[[rory]][1],colnames(movingbehavior))]!="stationary",1,0)
    leavingthismouse<-as.data.frame(leavingthismouse,ncol=1)
    colnames(leavingthismouse)<-  gsub(":",".lvs.",gsub("-angle","",OrderedMicePairs))[rory]

    if(rory<2){
      leavingbehavior<-leavingthismouse
    } else {
      leavingbehavior<-cbind(leavingbehavior,leavingthismouse)
    }
  }
  
  oknow<-cbind(RunMouse,movingbehavior,approachingbehavior,leavingbehavior)
  
  return(oknow)
}





# Mouse2MouseTrajectories -------------------------------------------------
# Function that takes two columns, per individual, of x/y coordinates
# (e.g. 4 individuals, n.inds=4, should correspond to 8 columns (xyvalues.allpairs))
# Returns pairwise distances, orientations at each time point 
# (i.e. whether a mouse is moving towards another) and velocities
Mouse2Mouse<-function(xyvalues.allpairs,pairwisedistances,
                      n.inds=4,shrinksize=.754){ 
  #xyvalues.allpairs = subset of full dataframe containing original all x,y coordinates for all
  # n.inds = number of individual mice
 
  
  ncomparisons<-2*(n.inds-1)
  flag<-1
  colorlist<-c("gray","blue","red","green")

    for(d in 1:n.inds){ #Loops through all individuals, to compare trajectory relative to others location
     
      
      
      a1<-(d-1)*2+1
      b1<-a1+1
      Z1=xyvalues.allpairs[,a1]+1i*xyvalues.allpairs[,b1] #Pulls locations for first individual
      # RealZ <-zoo(Z1)
      # #Location Vector 
      # Z1<-(RealZ)
      
      # step vectors
      dZ1 <- (diff(Z1))
      # orientation of each step
      distanctrav<-Mod(dZ1)
      Z1.Phi <- Arg(dZ1)
      Compass.D.M1<-((Z1.Phi * 180)/pi)
      Compass.D.M1<-ifelse(Compass.D.M1>0,Compass.D.M1,(360+(Compass.D.M1)))
        
      
      circle.Compass.D.M1<-circular(Compass.D.M1, type = "angles",
               units = "degrees",
               template = "geographics",#c("none", "geographics", "clock12", "clock24"),
               #modulo = c("asis", "2pi", "pi"),
               zero = 0) #rotation = c("counter", "clock"), names)
      
     nom1<-gsub(".x","",colnames(xyvalues.allpairs)[a1])
     # rose.diag(circle.Compass.D.M1, bins = 16, col = colorlist[d], 
     #           prop = 2, shrink=shrinksize,
     #           main = nom1,rotation="clock")

      
        for(e in 1:n.inds){  #Loops through all others
          if(d!=e){
          
          a2<-(e-1)*2+1
          b2<-a2+1
          
          
          comparisonheader<-paste(gsub(".x","",colnames(xyvalues.allpairs)[a1]),
                                  gsub(".x","",colnames(xyvalues.allpairs)[a2]),sep=':')
          
          
          Z2=xyvalues.allpairs[,a2]+1i*xyvalues.allpairs[,b2] #Pulls locations for first individual
 
          #interleave movement and difference vectors
          interleaveddata<-matrix(rbind(t(Z1), t(Z2)), ncol=1, byrow=TRUE)
          
          # step vectors
          dZ.1.2 <- diff(interleaveddata)#Calculates distance between 1 & 2
          
          # orientation of each step
          Z1Z2.Phi <- Arg(dZ.1.2)#orientation of vector connecting mouse 1 and mouse 2
          Z1Z2.Phi <- Z1Z2.Phi[seq(1,nrow(Z1Z2.Phi),2),]
          Z1Z2.Phi <- Z1Z2.Phi[c(1:(length(Z1Z2.Phi)-1))]
          
          Compass.D.M1toM2<-((Z1Z2.Phi * 180)/pi)#Calculates orientation between M1 and M2, in angles
          Compass.D.M1toM2<-ifelse(Compass.D.M1toM2>0,Compass.D.M1toM2,(360+(Compass.D.M1toM2)))
          
          
          # AngularDifference<-atan2(sin(Z1Z2.Phi-Z1.Phi), cos(Z1Z2.Phi-Z1.Phi))
          # DiffCompass.D.M1toM2<-(90 - ((AngularDifference) * 180)/pi)#Calculates orientation between M1 and M2, in angles
          # DiffCompass.D.M1toM2<-ifelse(DiffCompass.D.M1toM2>0,DiffCompass.D.M1toM2,(360+(DiffCompass.D.M1toM2)))#makes angles go 0-360
          # 
          
          #Two-step function to calculate difference between two angles (which might span 0), and for
          # which negative values are returned for one direction of turn, and positive for another
          a = Compass.D.M1toM2 - Compass.D.M1
          a = ifelse(a>180,a-360,
                     ifelse(a<(-180),a+360,a))
          


          #up, down, right, left,upright,upleft,downright,downleft
          # (defines pairwise angle for plotting quadrant to quadrant angles with relevant orientation)
          anglecomparisons<-c(pi/2,1.5*pi,0,pi,
                              pi/4,.75*pi,pi+.75*pi,pi+pi/4)
          
          #Depending on which quadrants are being compared, this nested ifelse will assign zervalues corresponding to the 
          # different elements in the anglecomparisons vector
          zervalue<-ifelse(comparisonheader=="A3:A1" || comparisonheader=="A4:A2",anglecomparisons[1],
                           ifelse(comparisonheader=="A1:A3" || comparisonheader=="A2:A4",anglecomparisons[2],
                                  ifelse(comparisonheader=="A1:A2" || comparisonheader=="A3:A4",anglecomparisons[3],
                                         ifelse(comparisonheader=="A4:A3" || comparisonheader=="A2:A1", anglecomparisons[4],
                                                ifelse(comparisonheader=="A3:A2",anglecomparisons[5],
                                                       ifelse(comparisonheader=="A4:A1",anglecomparisons[6],
                                                              ifelse(comparisonheader=="A1:A4",anglecomparisons[7],anglecomparisons[8])))))))
          
          #Turns angle data into a 'circular' object, for circular plotting
          circle.Compass.M1.M2<-circular(a, type = "angles",
                                        units = "degrees",
                                        template = "none",#c("none", "geographics", "clock12", "clock24"),
                                        #modulo = c("asis", "2pi", "pi"),
                                        zero = zervalue) #rotation = c("counter", "clock"), names)
          
          
          
          rose.diag(circle.Compass.M1.M2, bins = 16, col = colorlist[e], prop = 2, shrink=shrinksize,
                    main =comparisonheader,col.lab=colorlist[d] ,rotation="clock")
          
          
          
          #From Almelda et al. 2010, Indices of movement behaviour: conceptual background, effects of scale and location errors.
          # "The Straightness or linearity index, ST, (BATSCHELET 1981), 
          # is simply the net displacement distance (the Euclidian distance 
          # between the start and the final point), divided by the total length of the movement."
          
          #BATSCHELET, E. 1981. Circular Statistics in Biology. London, Academic Press.
          
          
          if(flag>1){
            allcomparisonangles<-cbind(allcomparisonangles,a)
            colnames(allcomparisonangles)[ncol(allcomparisonangles)]<-paste(comparisonheader,"angle",sep='-')
          } else {
            allcomparisonangles<-as.data.frame(a,ncol=1)
            colnames(allcomparisonangles)<-paste(comparisonheader,"angle",sep='-')
          }
          
          
          
          
          #############################
          # TotDistance<-sum(S) #Sum of step lengths
          # AvgSpeed<-mean(V)   #Mean of velocities
          # MaxSpeed<-max(V)    #Maximum velocity
          # StraightnessIndex<-(StrtStop/TotDistance) #Ranges from 0-1 (with 1 being a straight line)
          # AvgOrientation<-mean(Phi)
          # AvgTurn<-mean(Theta)
          
          

          flag<-flag+1
          } else {
            rose.diag(circle.Compass.D.M1, bins = 16, col = colorlist[d], 
                      prop = 2, shrink=shrinksize,
                      main = nom1,rotation="clock")
          }
        }  
      
      
     if(d>1){
       maintraveldis<-cbind(maintraveldis,distanctrav)
       colnames(maintraveldis)[ncol(maintraveldis)]<-paste(nom1,"step",sep='')
     } else {
       maintraveldis<-as.data.frame(distanctrav,ncol=1)
       colnames(maintraveldis)<-paste(nom1,"step",sep='')
     }
     
      if(d>1){
        maintravelangles<-cbind(maintravelangles,Compass.D.M1)
        colnames(maintravelangles)[ncol(maintravelangles)]<-paste(nom1,"angles",sep='-')
      } else {
        maintravelangles<-as.data.frame(Compass.D.M1,ncol=1)
        colnames(maintravelangles)<-paste(nom1,"angles",sep='-')
      }
    }
  
    allcomparisons<-cbind(maintraveldis,maintravelangles,allcomparisonangles)
  
    
  #Add on data.frame construction, making new columns corresponding 
  # to change (Delta) in pairwise distances (to be used when using movement data for behavioral categorization)
    G<-apply(pairwisedistances,2,function(x) c(0,diff(x)))
    colnames(G)<-paste("Delta.",colnames(pairwisedistances),sep='')
    G<-G[-1,]
    
    allcomparisons<-cbind(allcomparisons,G)
    

  
  
  
  return(allcomparisons)
}






# pairwise.Distances ------------------------------------------------------
#Function to calculate pairwise distances from CombinedInfo dataframes
# *dataframeCombinedInfo* should be a df (nrows corresponding to frames), with labelled x and y coordinates 
# (e.g. "x.a1"    "y.a1"    "x.a2"    "y.a2"    "x.a3"    "y.a3"    "x.a4"    "y.a4" )
# *inds* should correspond to the number of individuals you are comparing


pairwise.Distances<-function(dataframeCombinedInfo,inds=4,pix.cm=4.2924){
  if(inds>1){
  individual.names<-gsub(".x","",colnames(dataframeCombinedInfo)[seq(2,inds*2,2)])#pulls subset of column names, corresponding to the # of inds, and uses gsub to substitute 'blank' for x.
  ncomparisons<-2*(inds-1) #Calculates number of unique dyadic comparisons based on the n of individuals
  for(i in 1:(inds-1)){
    first<-dataframeCombinedInfo[,grepl(individual.names[i],colnames(dataframeCombinedInfo))]
    firstID<-individual.names[i]
    for(j in (i+1):inds){
      second<-dataframeCombinedInfo[,grepl(individual.names[j],colnames(dataframeCombinedInfo))]
      distances<-sqrt((first[,1]-second[,1])^2+(first[,2]-second[,2])^2)
      distances<-distances/pix.cm
      secondID<-individual.names[j]
      comparisonname<-paste(firstID,secondID,"dist",sep='')
      dataframeCombinedInfo<-cbind(dataframeCombinedInfo,distances)
      colnames(dataframeCombinedInfo)[ncol(dataframeCombinedInfo)]<-comparisonname
    }
  }
  return(dataframeCombinedInfo) #returns original dataframe with additional columns corresponding to dyadic distances at each timepoint
  }
}



# Associate.Identifier ----------------------------------------------------
#Using the *dataframeCombinedInfo*, which should contain pairwise distance columns with "dist" in the column names
# cheks for associations based on whether the distances are smaller than 'AssociatingDistance' threshold

Associate.Identifier<-function(dataframeCombinedInfo,AssociatingDistance){
  associationspossible<-colnames(dataframeCombinedInfo[,grepl("dist",colnames(dataframeCombinedInfo))])
  for(r in 1:length(associationspossible)){
    ass<-associationspossible[r]
    ass.lab<-gsub("dist","",ass)
    xyz<-ifelse(dataframeCombinedInfo[,ass]<AssociatingDistance,1,0)
    dataframeCombinedInfo<-cbind(dataframeCombinedInfo,xyz)
    colnames(dataframeCombinedInfo)[ncol(dataframeCombinedInfo)]<-ass.lab
  }
  #CombinedInfo$totalSocial<-apply(CombinedInfo[, c(16:21)], 1, function(x) toString(na.omit(x)))
  return(dataframeCombinedInfo)
}



####
# Takes HugeMouse dataframe (with colnames id'd below),
# where 'ds' represents dailysecond, where each ds is repeated 10x
# and creates a downsampled version where each set of measures collected over each
# of the 10 frames/second are used to generate metrics at 1 Htz
#
# For categorical movement states (running, walking, stationary), takes the most freq state
# For velocity (cm/s), takes average over 10 frames/sec
# For binary (0/1) variables, if it happens at all, take it!
###################################################################
# c("ds", "frame", "a1.x", "a1.y", "a2.x", "a2.y", "a3.x", "a3.y", 
#   "a4.x", "a4.y", "a1a2dist", "a1a3dist", "a1a4dist", "a2a3dist", 
#   "a2a4dist", "a3a4dist", "a1a2", "a1a3", "a1a4", "a2a3", "a2a4", 
#   "a3a4", "a1step", "a2step", "a3step", "a4step", "a1-angles", 
#   "a2-angles", "a3-angles", "a4-angles", "a1:a2-angle", "a1:a3-angle", 
#   "a1:a4-angle", "a2:a1-angle", "a2:a3-angle", "a2:a4-angle", "a3:a1-angle", 
#   "a3:a2-angle", "a3:a4-angle", "a4:a1-angle", "a4:a2-angle", "a4:a3-angle", 
#   "Delta.a1a2dist", "Delta.a1a3dist", "Delta.a1a4dist", "Delta.a2a3dist", 
#   "Delta.a2a4dist", "Delta.a3a4dist", "a1cm/s", "a2cm/s", "a3cm/s", 
#   "a4cm/s", "a1movement", "a2movement", "a3movement", "a4movement", 
#   "a1.app.a2", "a1.app.a3", "a1.app.a4", "a2.app.a1", "a2.app.a3", 
#   "a2.app.a4", "a3.app.a1", "a3.app.a2", "a3.app.a4", "a4.app.a1", 
#   "a4.app.a2", "a4.app.a3", "a1.lvs.a2", "a1.lvs.a3", "a1.lvs.a4", 
#   "a2.lvs.a1", "a2.lvs.a3", "a2.lvs.a4", "a3.lvs.a1", "a3.lvs.a2", 
#   "a3.lvs.a4", "a4.lvs.a1", "a4.lvs.a2", "a4.lvs.a3")
HugeMouseDownSampler<-function(HugeMouse){
  
  numericcolumnstoavg<-colnames(HugeMouse)[c(grep("cm/s",colnames(HugeMouse)))]
  numericcolumnstomax<-colnames(HugeMouse)[c(which(colnames(HugeMouse) %in% c("a1a2","a1a3","a1a4","a2a3","a2a4","a3a4")),
                         grep(".app.",colnames(HugeMouse)),grep("lvs",colnames(HugeMouse)))]
  
  movcols<-colnames(HugeMouse)[grep("movement",colnames(HugeMouse))]
  

  
  for(nca in 1:length(numericcolumnstoavg)){
    avgspd<-HugeMouse %>%
      group_by(ds) %>%
        summarise_at(numericcolumnstoavg[nca], mean, na.rm = TRUE)
    #############################
      if(nca==1){
        howmove<-avgspd
      } else {
        howmove<-merge(howmove,avgspd,by="ds")
      }
  }
  
  for(ncm in 1:length(numericcolumnstomax)){
    maxVv<-HugeMouse %>%
      group_by(ds) %>%
      summarise_at(numericcolumnstomax[ncm], max, na.rm = TRUE)
    #############################
    howmove<-merge(howmove,maxVv,by="ds")
    
  }
  
  for(mc in 1:4){
    howmove[,movcols[mc]]<-ifelse(howmove[,numericcolumnstoavg[mc]]>10,"running",
                                  ifelse(howmove[,numericcolumnstoavg[mc]]>0.5,"walking","stationary"))
  }
  
  howmove<-howmove[,c(1,36:39,2:35)]
  
  howmove<-myFun(howmove)
  
  howmove[,c(2:5)]<-lapply(howmove[,c(2:5)],function(x){factor(x)})
  
  return(howmove)
}



#-----------------------------------------------------------------------------------
# Function that takes a vector of movement 'states' ('running','walking','stationary')
# and returns the single value that is most represented 
# (or, if running is tied for most represented, takes 'running' as the state), else it looks at walking
# and last, if neither running nor walking is top (or tied for top), then the state is
# 'stationary'
run.walk.station<-function(t){
  t2<-data.frame(table(t))
  rn<-t2[which(t2$t=="running"),'Freq']
  
  #if most prevalent beh state is running, or tied for most prevalent, then it's running
  if(max(t2$Freq)==rn){
    bstate<-'running'
  } else {
    wk<-t2[which(t2$t=="walking"),'Freq']
    #if most prevalent beh state is walking, then it's walking
    if(max(t2$Freq)==wk){
      bstate<-'walking'
    } else {
      bstate<-'stationary'
    }
  }
  return(bstate)
}

#############
myFun <- function(data) {
  ListCols <- sapply(data, is.list)
  cbind(data[!ListCols], t(apply(data[ListCols], 1, unlist)))
}#FUNCTION TO UNLIST VARIALBES IN DF
########################################################

# add.alpha ---------------------------------------------------------------
## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}



# characteriseTrajectory --------------------------------------------------

# Define a function which calculates some statistics
# of interest for a single trajectory
characteriseTrajectory <- function(trj) {
  # Measures of speed
  derivs <- TrajDerivatives(trj)
  mean_speed <- mean(derivs$speed)
  sd_speed <- sd(derivs$speed)
  
  # Measures of straightness
  sinuosity <- TrajSinuosity(trj)
  resampled <- TrajRediscretize(trj, .001)
  Emax <- TrajEmax(resampled)
  
  # Periodicity
  corr <- TrajDirectionAutocorrelations(resampled, 60)
  first_min <- TrajDAFindFirstMinimum(corr)
  
  # Return a list with all of the statistics for this trajectory
  list(mean_speed = mean_speed,
       sd_speed = sd_speed,
       sinuosity = sinuosity,
       Emax = Emax,
       min_deltaS = first_min[1],
       min_C = first_min[2]
  )
}




# chart.Correlation.RUSTY --------------------------------------------------
#Modified from PerformanceAnalytics::chart.Correlation
chart.Correlation.RUSTY<-function (R, histogram = TRUE, method = c("pearson", "kendall", 
                                                                   "spearman"), dotcolors, ppp, logger, ...) 
{
  x = checkData(R, method = "matrix")
  if (missing(method)) 
    method = method[1]
  cormeth <- method
  panel.cor <- function(x, y, digits = 2, prefix = "", 
                        use = "pairwise.complete.obs", method = cormeth, 
                        cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 0.8/strwidth(txt)
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
                                                                              "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
    text(0.8, 0.8, Signif, cex = cex, col = 2)
  }
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  hist.panel = function(x, ... = NULL) {
    par(new = TRUE)
    hist(x, col = "light gray", probability = TRUE, 
         axes = FALSE, main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    rug(x)
  }
  
  if(missing(logger))
    logger=""
  
  if (histogram) 
    pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
          diag.panel = hist.panel,col=dotcolors,pch=ppp,log=logger)
  else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor,col=dotcolors,pch=ppp,log=logger)
}
