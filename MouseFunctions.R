












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
