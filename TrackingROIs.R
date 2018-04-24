




setwd("C:/Users/Rusty/Amazon Drive/MICE")
source("MouseFunctions.R", chdir = F)


library(MASS)


setwd(dir<-"C:/Users/Rusty/Amazon Drive/MICE/IR")#sets working directory (setwd), and new variable (dir) corresponding to parent directory


directories<-list.dirs(dir,recursive = TRUE)#list.dirs lists directories 
directories<-directories[-1]#removes "self" directory




# CreateCompositeList -----------------------------------------------------
#Loops through each folder in 'directories'
# pulls files corresponding to tracking informaiton (location.names)
# and ROI information (Roifiles)
# then names the information from these files and amalgamates into a data.frame called CombinedInfo
# then puts that dataframe into the list AllFoldersList at position zz

AllFoldersList<-list() #Creates empty list which will get filled iteratively with loop
pix.cm<-4.2924
AssociatingDistance.cm<-10 # Define cm distance that corresponding to an 'association'
flag<-1 #sets flag, which will only increase when looping through folders containing the "correct" kind of data
for (zz in 1:length(directories)){
  FolderInfo1<-directories[zz]#pull full 
  FolderInfo<-strsplit(FolderInfo1,"/")[[1]][length(strsplit(FolderInfo1,"/")[[1]])]
  FolderInfo<-paste(strsplit(FolderInfo1,"/")[[1]][(length(strsplit(FolderInfo1,"/")[[1]])-1)],FolderInfo,sep='.')
  location.names<-list.files(directories[zz],full.names=TRUE,pattern="fixed.csv")#change pattern to only pull location csvs
  
  if(length(location.names)>0){ #Only if in the right kind of folder, i.e. containing ...fixed.csv files, run the rest, otherwise, skip
  
    Roifiles<-list.files(directories[zz],full.names = TRUE,pattern="-region")#refine to pull roi csvs
    Roifilenames<-list.files(directories[zz],full.names = FALSE,pattern="-region")#refine to pull roi csvs
    inds<-length(location.names)
    
    for(Caleb in 1:inds){
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
    #CombinedInfo<-round(CombinedInfo)
    CombinedInfo<-pairwise.Distances(CombinedInfo,inds)#Custom function located in MouseFunctions.R
    ncomparisons<-2*(inds-1) #Calculates number of unique dyadic comparisons based on the n of individuals
    dister<-(inds*2+2)
    CombinedInfo[,c(dister:(dister+ncomparisons-1))]<-(CombinedInfo[,c(dister:(dister+ncomparisons-1))])/(pix.cm)
    CombinedInfo<-Associate.Identifier(CombinedInfo,AssociatingDistance.cm)#Custom function located in MouseFunctions.R
    
    for(Caitlin in 1:length(Roifiles)){
      ROI<-read.csv(Roifiles[Caitlin])
      roiname<-Roifilenames[Caitlin]
      colnames(ROI)[2]<-roiname
      CombinedInfo<-merge(CombinedInfo,ROI, by="position")
    }
    
    
    AllFoldersList[[flag]]<-CombinedInfo #puts CombinedInfo dataframe into AllFoldersList at position zz
    names(AllFoldersList)[[flag]]<-FolderInfo #applies name of folder to list element
    flag<-flag+1
  }
}  











# TrialSpecificManipulations ----------------------------------------------

CombinedInfo<-AllFoldersList[[1]]#pull dataframe corresponding to particular day/trial
#CombinedInfo<-CombinedInfo[c(1:(10*floor(nrow(CombinedInfo)/10))),] #makes nrow divisible by 10 for subsampling purposes
#CombinedInfo.sub<-cbind(0.1*CombinedInfo[seq(1,length(CombinedInfo$position),10),1],
#                        apply(CombinedInfo[,c(2:21)],2,function(x) colMeans(matrix(x, nrow=10)))) #
# colnames(CombinedInfo.sub)[1]<-"position"
# CombinedInfo.sub<-as.data.frame(CombinedInfo.sub)
# summary(CombinedInfo.sub)



BASE<-matrix(summary(CombinedInfo),nrow=7)
colnames(BASE)<-colnames(CombinedInfo)



library(plyr)
check<-apply(CombinedInfo[,c(22:25)],2,function(x) plyr::count(x))
for(g in 1:length(check)){
  j<-check[[g]]
  colnames(j)[2]<-names(check)[g]
  if(g>1){
    boop<-cbind(boop,j)
  } else {
    boop<-j
  }
}



library(trajr)

A1trj <- TrajFromCoords(CombinedInfo[,c(2,3,1)],fps=10);
A2trj <- TrajFromCoords(CombinedInfo[,c(4,5,1)],fps=10)
A3trj <- TrajFromCoords(CombinedInfo[,c(6,7,1)],fps=10)
A4trj <- TrajFromCoords(CombinedInfo[,c(8,9,1)],fps=10)

All.arenas<-list(A1trj,A2trj,A3trj,A4trj)

for(tjmax in 1:length(All.arenas)){
  sp.Traj<-All.arenas[[tjmax]]
  print(paste("Trajectory length ",TrajLength(sp.Traj)/100/1000/pix.cm))
  sp.Traj$displaceDist<-c(0,Mod(diff(sp.Traj$polar[0:nrow(sp.Traj)])))
  sorted.displace<-sp.Traj[order(-sp.Traj$displaceDist),]
}

ch<-c(NA,TrajAngles(sp.Traj,compass.direction = 0))




# Calculate all stats for trajectories in the list
# which was built in the previous example
##stats <- TrajsMergeStats(All.arenas, characteriseTrajectory) #currently this command/function is not optimized and takes too long (>5 hours)
#print(stats)
beepr::beep(4)






# Calculate speed and acceleration
derivs <- TrajDerivatives(A1trj)










#pdf("A1movement_long.pdf",width= 8, height= 2,family="NimbusRom")
png("A1movement_long.png",width= 6400, height= 800,family="NimbusRom")

par(mfrow=c(1,1))# rows, columns
par(mar=c(6,5,0,6)) #Margines of each plot (bottom, left, top, right)
#par(oma=c(0,0,10,0)) #Outer Margins of "entire" plot
# Plot acceleration and speed
plot(derivs$acceleration ~ derivs$accelerationTimes, type = 'l', col = 'red', 
     yaxt = 'n',
     xlab = 'Time (s)',
     ylab = expression(paste('Acceleration (', m/s^2, ')')))
axis(side = 2, col = "red")
lines(derivs$speed ~ derivs$speedTimes, col = 'blue')
axis(side = 4, col = "blue")
mtext('Speed (m/s)', side = 4, line = 3)
abline(h = 0, col = 'lightGrey')
dev.off()


TrajSpeedIntervals
TrajLength(A1trj)





#PLOTTING

#Total movement
plot(1,type="n", xlab="", ylab="",ylim=c(0,1080),xlim=c(0,1400))
lines(y0.A1~x0.A1,data=CombinedInfo,ylim=c(0,1080),xlim=c(0,1400),col=add.alpha("gray",0.8))
lines(y0.A2~x0.A2,data=CombinedInfo,ylim=c(0,1080),xlim=c(0,1400),col=add.alpha("blue",0.8))
lines(y0.A3~x0.A3,data=CombinedInfo,ylim=c(0,1080),xlim=c(0,1400),col=add.alpha("red",0.8))
lines(y0.A4~x0.A4,data=CombinedInfo,ylim=c(0,1080),xlim=c(0,1400),col=add.alpha("green",0.8))


#All mice, needed for heatmaps or 3d plots
allheat<-kde2d(c(CombinedInfo$x0.A1,CombinedInfo$x0.A2,CombinedInfo$x0.A3,CombinedInfo$x0.A4),
              c(CombinedInfo$y0.A1,CombinedInfo$y0.A2,CombinedInfo$y0.A3,CombinedInfo$y0.A4), h=75, n=200)
filled.contour(allheat,color.palette=colorRampPalette(c('purple','blue','yellow','red','darkred')),ylim=c(0,1080),xlim=c(0,1400))

colfunc<-colorRampPalette(rev(c("red","yellow","springgreen","royalblue","black")))
image(allheat$x, allheat$y, allheat$z, col = colfunc(50), axes = TRUE,ylab="Vertical coordinates",
      xlab="Horizontal coordinates",main="Activity")

require(rgl)
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(allheat$z, nbcol)
persp3d(x=allheat$x,y=allheat$y,
        z = allheat$z, theta = 120,col=color[zcol],
        zlab="Time spent",
        xlab="Horizontal coordinates",ylab="Vertical coordinates",main="Activity patterns")

#######################

colfunc.a1<-colorRampPalette(rev(c("darkgoldenrod1","black","white")))
colfunc.a2<-colorRampPalette(rev(c("purple","blue","cyan")))
colfunc.a3<-colorRampPalette(rev(c("mediumorchid1","orange","yellow")))
colfunc.a4<-colorRampPalette(rev(c("dodgerblue","darkgreen","mediumspringgreen")))


a1heat<-kde2d(CombinedInfo$x0.A1, CombinedInfo$y0.A1, h=75, n=150)
a2heat<-kde2d(CombinedInfo$x0.A2, CombinedInfo$y0.A2, h=75, n=150)
a3heat<-kde2d(CombinedInfo$x0.A3, CombinedInfo$y0.A3, h=75, n=150)
a4heat<-kde2d(CombinedInfo$x0.A4, CombinedInfo$y0.A4, h=75, n=150)


image(a1heat,ylim=c(0,1080),xlim=c(0,1400),col = colfunc(50))
par(new = TRUE)
image(a2heat,ylim=c(0,1080),xlim=c(0,1400),col = colfunc(50))
par(new = TRUE)
image(a3heat,ylim=c(0,1080),xlim=c(0,1400),col = colfunc(50))
par(new = TRUE)
image(a4heat,ylim=c(0,1080),xlim=c(0,1400),col = colfunc(50))
lines(y0.A1~x0.A1,data=CombinedInfo,ylim=c(0,1080),xlim=c(0,1400),col=add.alpha("gray",0.2))
lines(y0.A2~x0.A2,data=CombinedInfo,ylim=c(0,1080),xlim=c(0,1400),col=add.alpha("blue",0.2))
lines(y0.A3~x0.A3,data=CombinedInfo,ylim=c(0,1080),xlim=c(0,1400),col=add.alpha("red",0.2))
lines(y0.A4~x0.A4,data=CombinedInfo,ylim=c(0,1080),xlim=c(0,1400),col=add.alpha("darkgreen",0.2))




zcol  = cut(a1heat$z, nbcol)
persp3d(x=a1heat$x,y=a1heat$y,
        z = a1heat$z, theta = 120,col=color[zcol],ylim=c(0,1080),xlim=c(0,1400),
        zlab="Time spent",
        xlab="Horizontal coordinates",ylab="Vertical coordinates",main="Activity patterns")
zcol  = cut(a2heat$z, nbcol)

persp3d(x=a2heat$x,y=a2heat$y,
        z = a2heat$z, theta = 120,col=color[zcol],ylim=c(0,1080),xlim=c(0,1400),
        zlab="Time spent",
        xlab="Horizontal coordinates",ylab="Vertical coordinates",main="Activity patterns", add=TRUE) 
zcol  = cut(a3heat$z, nbcol)

persp3d(x=a3heat$x,y=a3heat$y,
        z = a3heat$z, theta = 120,col=color[zcol],ylim=c(0,1080),xlim=c(0,1400),
        zlab="Time spent",
        xlab="Horizontal coordinates",ylab="Vertical coordinates",main="Activity patterns", add=TRUE) 
zcol  = cut(a4heat$z, nbcol)

persp3d(x=a4heat$x,y=a4heat$y,
        z = a4heat$z, theta = 120,col=color[zcol],ylim=c(0,1080),xlim=c(0,1400),
        zlab="Time spent",
        xlab="Horizontal coordinates",ylab="Vertical coordinates",main="Activity patterns", add=TRUE) 
