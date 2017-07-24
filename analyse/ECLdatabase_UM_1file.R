##This is code to take the output from the low tracking software 
##and turn it into a .csv file of ECL track data.

##Note that this script performs only the initial filtering 
##to restrict ECLs based on location, track length and at least one closed fix
##To obtain a better database of ECLs with impacts you can filter for 
##a fix curvature of >= 0.25 & the fix being in the ECL region

##This version of the code requires both the annual version of the data
##as well as that run for the previous summer, to collate the tracks for events
##over January 1.

##Note also that R has issues with large datasets - 
##I tend to separate into sets of 20 years or fewer. Pretty quick. 

##yearS and yearE are start & end years. Must have DJF data for yearS-1
##output is a string to help label the fix & event files (so use "name")

##You will need to set the location to where the track files are stored.
#setwd('~/Charles_bom/output')

ECLdatabase<-function(fname,year)
{

  read.table(paste("tracks_",fname,".dat",sep=""), sep="",skip=0)->data
  date1=data[,3]
  dateB=date1 %% 10000
  dateA=floor(date1/10000)
  dlist=unique(dateA)
  data[,3]=year*10000+dateB
  fixes<-data

  ##Identify if within the ECL region  
  fixes[,14]<-0
  I<-which(fixes[,6]>=149 & fixes[,6]<=161 & fixes[,7]<(-37) & fixes[,7]>=-41)
  fixes[I,14]<-1
  I<-which(fixes[,6]>=(149+(37+fixes[,7])/2) & fixes[,6]<=161 & fixes[,7]<(-31) & fixes[,7]>=-37)
  fixes[I,14]<-1
  I<-which(fixes[,6]>=152 & fixes[,6]<=161 & fixes[,7]<=(-24) & fixes[,7]>=-31)
  fixes[I,14]<-1
    
  ##Collates information for events
  x<-rle(fixes[,1])
  events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=8))
  
  ##Single-fix events
  
  I<-which(events[,2]==1)
  if(length(I)>0)
  {
    y<-match(events[I,1],fixes[,1])
    events=events[-I,]
    fixes=fixes[-y,]
  }
  
  for(i in 1:length(events[,1])) 
  {
    I<-which(fixes[,1]==events[i,1])
    events[i,3]=min(fixes[I,3]) ##Date1
    events[i,4]=max(fixes[I,3]) ##Date2
    events[i,5]=min(fixes[I,8]) ##Min central pressure
    events[i,6]=max(fixes[I,9]) ##Max curvature
    events[i,7]=min(fixes[I,2]) ##First fix # - for cross-year events
    events[i,8]=max(fixes[I,2]) ##Last fix #
    if(sum(fixes[I,14])>0) events[i,9]=1 ##Ever in ECL region
    events[i,10]<-min(fixes[I,5] %% 10) ##Ever closed  
  }
  
  ## Take only the events where at least one fix is closed & in the ECL region
  ##And duration is at least 2 consecutive fixes
  
  I<-which(events[,9]>0 & events[,2]>1 & events[,10]==0)
  events2<-events[I,]
  include<-match(fixes[,1],events2[,1])
  J<-which(is.na(include)==0)
  fixes2=fixes[J,]

  ##Fix numbering and save final data
  
  eventnum=events2[,1]
  
  for(i in 1:length(eventnum))
  {
    I<-which(fixes2[,1]==eventnum[i])
    fixes2[I,1]=i
    events2[i,1]=i
  }

  ##Convert to dataframe
  
  Events<-data.frame(events2[,1:6])
  names(Events)=c('ID','Length','Date1','Date2','MSLP','CV')
  
  Fixes<-data.frame(fixes2[,1:9],fixes2[,14],fixes2[,11:13])
  names(Fixes)=c('ID','Fix','Date','Time','Open','Lon','Lat','MSLP','CV','Location','Radius','U movement','V movement')
  
  ##Add extra columns - count,mslp,cv,rad in loc
  Events$Rad2<-Events$CV2<-Events$MSLP2<-Events$Length2<-rep(0,length(Events$ID))

  for(i in 1:length(Events[,1]))
  {
    I<-which(Fixes$ID==Events[i,1] & Fixes$Location==1)
    Events$Length2[i]=length(I)
    Events$MSLP2[i]=min(Fixes$MSLP[I])
    Events$CV2[i]=max(Fixes$CV[I])
    Events$Rad2[i]=mean(Fixes$Radius[I])
  }

  outF=paste('ECLfixes_',fname,'.csv',sep="")
  outE=paste('ECLevents_',fname,'.csv',sep="")
  
  write.csv(Fixes,outF)
  write.csv(Events,outE)
}