


library(gdata)

# Mouse2MouseTrajectories -------------------------------------------------

Mouse2Mouse<-function(xyvalues.allpairs,n.inds=4,interval.of.frames=0.1){ 
  #xyvalues.allpairs = subset of full dataframe containing original all x,y coordinates for all
  # n.inds = number of individual mice
 # interval.of.frames used for calculating velocity towards one another (0.1 = 10fps)
  
  ncomparisons<-2*(n.inds-1)
  flag<-1
  
    for(d in 1:n.inds){ #Loops through all individuals, to compare trajectory relative to others location
     
      a1<-(d-1)*2+1
      b1<-a1+1
      Z1=xyvalues.allpairs[,a1]+1i*xyvalues.allpairs[,b1] #Pulls locations for first individual
      RealZ <-zoo(Z1)
      #Location Vector 
      Z1=as.data.frame(RealZ)
      
        for(e in 1:inds){  #Loops through all others
  
          a2<-(e-1)*2+1
          b2<-a2+1
          
          
          comparisonheader<-paste(gsub("x0.","",colnames(xyvalues.allpairs)[a1]),
                                  gsub("x0.","",colnames(xyvalues.allpairs)[a2]),sep=':')
          
          
          Z2=xyvalues.allpairs[,a2]+1i*xyvalues.allpairs[,b2] #Pulls locations for first individual
          RealZ2 <-zoo(Z2)
          #Location Vector 
          Z2=as.data.frame(RealZ2)
          
  
          interleaveddata<-matrix(rbind(t(Z1), t(Z2)), ncol=1, byrow=TRUE)
          
          # step vectors
          dZ <- diff(interleaveddata)
          
          # orientation of each step
          Phi <- Arg(dZ)
          
          # turning angles
          Theta <- diff(Phi)
          
          #Step lengths
          S <- Mod(dZ)
          
          # time intervals
          #dT <- diff(Time1)
          
          # Magnitude of linear velocity between points
          V <- S/interval.of.frames
          
          #From Almelda et al. 2010, Indices of movement behaviour: conceptual background, effects of scale and location errors.
          # "The Straightness or linearity index, ST, (BATSCHELET 1981), 
          # is simply the net displacement distance (the Euclidian distance 
          # between the start and the final point), divided by the total length of the movement."
          
          #BATSCHELET, E. 1981. Circular Statistics in Biology. London, Academic Press.
          
          M2Mcomparison<-cbind(dZ,Phi,c(0,Theta),S,V)
          colnames(M2Mcomparison)<-paste(comparisonheader,c("StepVector","Orientation","TurningAngle","StepLength","Velocity"),sep='_')
          M2Mcomparison<-M2Mcomparison[seq(1,nrow(M2Mcomparison),2),]
          
          
          #############################
          # TotDistance<-sum(S) #Sum of step lengths
          # AvgSpeed<-mean(V)   #Mean of velocities
          # MaxSpeed<-max(V)    #Maximum velocity
          # StraightnessIndex<-(StrtStop/TotDistance) #Ranges from 0-1 (with 1 being a straight line)
          # AvgOrientation<-mean(Phi)
          # AvgTurn<-mean(Theta)
          
          
          if(flag>1){
            allcomparisons<-cbind(allcomparisons,M2Mcomparison)
          } else {
            allcomparisons<-M2Mcomparison
          }
          flag<-flag+1
        }  
      
    }
  return(allcomparisons)
}




# pairwise.Distances ------------------------------------------------------
#Function to calculate pairwise distances from CombinedInfo dataframes
# *dataframeCombinedInfo* should be a df (nrows corresponding to frames), with labelled x and y coordinates 
# (e.g. "x0.A1"    "y0.A1"    "x0.A2"    "y0.A2"    "x0.A3"    "y0.A3"    "x0.A4"    "y0.A4" )
# *inds* should correspond to the number of individuals you are comparing


pairwise.Distances<-function(dataframeCombinedInfo,inds=4){
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
