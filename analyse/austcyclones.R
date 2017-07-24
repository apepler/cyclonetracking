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

get_austevents<-function(yearS,yearE,output,type="high")
{

  ##Files with ECL track data for required years & preceding summers.  
  year<-seq(yearS,yearE)
  fname1=paste('tracks_',as.character(year),'.dat',sep="")

  ##collates all the annual tracks into one file
  ##And the extra ones at the start into a "Dec" file

  read.table(fname1[1], sep="",skip=0)->data
  yy=floor(data[,3]/10000)
  yy2=unique(yy)
  if(length(yy2)>1)
  {
    I=which(yy==yy2[1])
    data[I,3]=(data[I,3]%%10000) + (yearS-1)*10000
    fixesDEC=data[I,]
    I=which(yy==yy2[2])
    data[I,3]=(data[I,3]%%10000) + yearS*10000
    fixes<-data[I,]
  } else {
    data[,3]=(data[,3]%%10000) + yearS*10000
    fixes<-data[,]
    fixesDEC=matrix(NaN,1,length(data[1,]))
  }

  for(i in 2:length(year))
  {
    read.table(fname1[i], sep="",skip=0)->data
    yy=floor(data[,3]/10000)
    yy2=unique(yy)
    if(length(yy2)>1)
    {
     I=which(yy==yy2[1])
     data[I,3]=(data[I,3]%%10000) + (year[i]-1)*10000
     fixesDEC=rbind(fixesDEC,data[I,])

     I=which(yy==yy2[2])
     data[I,3]=(data[I,3]%%10000) + year[i]*10000
     fixes<-rbind(fixes,data[I,])
    } else {
     data[,3]=(data[,3]%%10000) + year[i]*10000
     fixes<-rbind(fixes,data[,])
    }
  }
  rm(data)

  ## Fix intensity if high?
  fixes[,9]=abs(fixes[,9])
  fixesDEC[,9]=abs(fixesDEC[,9])
  
  #Unique event numbers (starting with year)
  fixes[,1]<-fixes[,1]+floor(fixes[,3]/10000)*1000000
  
  ##Identify if within the aust region
  
  fixes[,14]<-0
  I<-which(fixes[,6]>=110 & fixes[,6]<=155 & fixes[,7]<(-10) & fixes[,7]>=-45)
  fixes[I,14]<-1
    
  fixesDEC[,1]<-fixesDEC[,1]+floor(fixesDEC[,3]/10000)*1000000
  fixesDEC[,14]<-0
  I<-which(fixesDEC[,6]>=110 & fixesDEC[,6]<=155 & fixesDEC[,7]<(-10) & fixesDEC[,7]>=-45)
  fixesDEC[I,14]<-1
  
  ##Collates information for events
  x<-rle(fixes[,1])
  events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=10))
  
  ##Single-fix events
  print("Single fix") 
  I<-which(events[,2]==1)
  y<-match(events[I,1],fixes[,1])
  events[I,3]<-events[I,4]<-fixes[y,3]
  events[I,5]<-fixes[y,8]
  events[I,6]<-fixes[y,9]
  events[I,7]<-events[I,8]<-fixes[y,2]
  events[I,9]<-fixes[y,14] 
  events[I,10]<-fixes[y,5] %% 10
  J=which(events[I,9]==1)
  events[I[J],11]=fixes[y[J],8]
  events[I[J],12]=fixes[y[J],9]
  
  ##Multi-fix events
  print("Multi-fix")
  J<-which(events[,2]>1)
  print(length(J)) 
  for(i in 1:length(J)) 
  {
    if((i %% 1000) == 0) print(i)
    I<-which(fixes[,1]==events[J[i],1])
    if(sum(fixes[I,14])>0)
    {
    events[J[i],3]=min(fixes[I,3]) ##Date1
    events[J[i],4]=max(fixes[I,3]) ##Date2
    if(type=="high") events[J[i],5]=max(fixes[I,8]) else events[J[i],5]=min(fixes[I,8]) ##Max central pressure
    events[J[i],6]=max(fixes[I,9]) ##Max curvature
    events[J[i],7]=min(fixes[I,2]) ##First fix # - for cross-year events
    events[J[i],8]=max(fixes[I,2]) ##Last fix #
    events[J[i],9]=1 ##Ever in ECL region
    events[J[i],10]<-min(fixes[I,5] %% 10) ##Ever closed
    events[J[i],11]=max(fixes[I,8]*fixes[I,14]) # All outside region = 0
    events[J[i],12]=max(fixes[I,9]*fixes[I,14]) # All outside region = 0

    }  
  }
  
  ## Try to identify cross-year events
  ## For each event that starts higher than 1, find the first row
  ## Then identify the row in the DJF file that has the same values for columns 2-6
  ## Next find a row in the fixes file that matches the row one higher in DJF
  ## Finally, apply the eventID of this event to all of the values for the second event
  print("Crossyear") 
  events2=events
  
  I<-which(events[,7]>1)
  ##Any case where the first fix is not fix #1
  I2=matrix(0,nrow=length(I),ncol=1)
  for(i in 1:length(I))
  {
    J<-which(fixes[,1]==events[I[i],1])
    K<-which(fixesDEC[,3]==fixes[J[1],3] & fixesDEC[,4]==fixes[J[1],4] & fixesDEC[,7]==fixes[J[1],7] & fixesDEC[,6]==fixes[J[1],6])
    L<-which(fixes[,3]==fixesDEC[K-1,3] & fixes[,4]==fixesDEC[K-1,4] & fixes[,7]==fixesDEC[K-1,7] & fixes[,6]==fixesDEC[K-1,6])
    if(length(L)==1)
    {
      ID<-fixes[L,1]
      fixes[J,1]<-ID
      I2[i]=1
      
      ##Make event1 represent the full event now
      M<-which(events[,1]==ID)
      events2[M,2]=events[M,2]+events[I[i],2]
      events2[M,4]=events[I[i],4]
      if(type=="high") events2[M,5]=max(events[I[i],5],events[M,5]) else events2[M,5]=min(events[I[i],5],events[M,5])
      events2[M,6]=max(events[I[i],6],events[M,6])
      events2[M,8]=events[I[i],8]
      events2[M,9]=max(events[I[i],9],events[M,9])
      events2[M,10]=min(events[I[i],10],events[M,10])
      if(type=="high") events2[M,11]=max(events[I[i],11],events[M,11]) else events2[M,11]=min(events[I[i],11],events[M,11])
      events2[M,12]=max(events[I[i],12],events[M,12])
 
    }
  }
  
  ##And remove the end-year events
  rm(fixesDEC)
  
  ## Take only the events where at least one fix is closed & in the region
  
  I<-which(events2[,9]>0 & events2[,10]==0)
  events3<-events2[I,]
  include<-match(fixes[,1],events3[,1])
  J<-which(is.na(include)==0)
  fixes2=fixes[J,]
  
  ##Fix numbering and save final data
  print("Tidying") 
  eventnum=events3[,1]
  
  for(i in 1:length(eventnum))
  {
    I<-which(fixes2[,1]==eventnum[i])
    fixes2[I,1]=i
    events3[i,1]=i
  }
  
  ##Convert to dataframe
  
  Events<-data.frame(events3[,c(1:6,11:12)])
  names(Events)=c('ID','Length','Date1','Date2','MSLP','CV',"Aust.MSLP","Aust.CV")
  
  Fixes<-data.frame(fixes2[,1:9],fixes2[,14])
  names(Fixes)=c('ID','Fix','Date','Time','Open','Lon','Lat','MSLP','CV','Location')
  
  outF=paste(output,'_fixes.csv',sep="")
  outE=paste(output,'_events.csv',sep="")
  
  write.csv(Fixes,outF)
  write.csv(Events,outE)
}

setwd('/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_highs_rad10cv0.075/')
get_austevents(1950,1970,"UM_highs_5070",type="high")
get_austevents(1970,1990,"UM_highs_7090",type="high")
get_austevents(1990,2016,"UM_highs_9016",type="high")

