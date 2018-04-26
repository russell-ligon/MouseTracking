


library(gdata)
library(zoo)
library(circular)




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
  colnames(movingbehavior)<-gsub("step","movement",colnames(movingbehavior))
  
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
    
    approachingthismouse<-ifelse(socialmovements[,rory]==1 & socialmovements[,c(grepl(oc3[[rory]][1],colnames(socialmovements))&
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
# Returns pairwise distances, orientations at each time point (i.e. whether a mouse is moving towards another) and velocities
Mouse2Mouse<-function(xyvalues.allpairs,pairwisedistances,n.inds=4,shrinksize=.75,pix.cm=4.2924){ 
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
      distanctrav<-Mod(dZ1)/pix.cm
      Z1.Phi <- Arg(dZ1)
      Compass.D.M1<-((Z1.Phi * 180)/pi)
      Compass.D.M1<-ifelse(Compass.D.M1>0,Compass.D.M1,(360+(Compass.D.M1)))
        
      
      circle.Compass.D.M1<-circular(Compass.D.M1, type = "angles",
               units = "degrees",
               template = "geographics",#c("none", "geographics", "clock12", "clock24"),
               #modulo = c("asis", "2pi", "pi"),
               zero = 0) #rotation = c("counter", "clock"), names)
      
     nom1<-gsub("x0.","",colnames(xyvalues.allpairs)[a1])
     rose.diag(circle.Compass.D.M1, bins = 16, col = colorlist[d], prop = 2, shrink=shrinksize,
                main = nom1,rotation="clock")

      
        for(e in 1:n.inds){  #Loops through all others
          if(d!=e){
          
          a2<-(e-1)*2+1
          b2<-a2+1
          
          
          comparisonheader<-paste(gsub("x0.","",colnames(xyvalues.allpairs)[a1]),
                                  gsub("x0.","",colnames(xyvalues.allpairs)[a2]),sep=':')
          
          
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
# (e.g. "x0.A1"    "y0.A1"    "x0.A2"    "y0.A2"    "x0.A3"    "y0.A3"    "x0.A4"    "y0.A4" )
# *inds* should correspond to the number of individuals you are comparing


pairwise.Distances<-function(dataframeCombinedInfo,inds=4){
  if(inds>1){
  individual.names<-gsub("x0.","",colnames(dataframeCombinedInfo)[seq(2,inds*2,2)])#pulls subset of column names, corresponding to the # of inds, and uses gsub to substitute 'blank' for x0.
  ncomparisons<-2*(inds-1) #Calculates number of unique dyadic comparisons based on the n of individuals
  for(i in 1:(inds-1)){
    first<-dataframeCombinedInfo[,grepl(individual.names[i],colnames(dataframeCombinedInfo))]
    firstID<-individual.names[i]
    for(j in (i+1):inds){
      second<-dataframeCombinedInfo[,grepl(individual.names[j],colnames(dataframeCombinedInfo))]
      distances<-sqrt((first[,1]-second[,1])^2+(first[,2]-second[,2])^2)
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
# cheks for associations based on whether the distances are smaller than 'AssociatingDistance.pixels' threshold

Associate.Identifier<-function(dataframeCombinedInfo,AssociatingDistance.pixels){
  associationspossible<-colnames(dataframeCombinedInfo[,grepl("dist",colnames(dataframeCombinedInfo))])
  for(r in 1:length(associationspossible)){
    ass<-associationspossible[r]
    ass.lab<-gsub("dist","",ass)
    xyz<-ifelse(dataframeCombinedInfo[,ass]<AssociatingDistance.pixels,1,0)
    dataframeCombinedInfo<-cbind(dataframeCombinedInfo,xyz)
    colnames(dataframeCombinedInfo)[ncol(dataframeCombinedInfo)]<-ass.lab
  }
  #CombinedInfo$totalSocial<-apply(CombinedInfo[, c(16:21)], 1, function(x) toString(na.omit(x)))
  return(dataframeCombinedInfo)
}



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
