## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

library(ncdf4)
extract_counts<-function(year1,year2,type="low",dir="/short/eg3/asp561/cts.dir/gcyc_out",thresh=0,dur=NA)
{
years=seq(year1,year2,1)
months=1:12

lat=seq(-89.5,89.5,1)
lon=seq(0.5,359.5,1)  ### Can always combine into bigger cells later
systems<-array(0,c(length(lon),length(lat),length(years),12))

for(y in 1:length(years))
{
print(years[y])
fname=paste(dir,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Year=floor(fixes$Date/10000)
fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=floor(fixes$Lat)
fixes$Lon2=floor(fixes$Lon)%%360
if(type=="high") fixes$CV=-fixes$CV

### Make table of events to combine with DJF for exclusion
 if(!is.na(dur))
  {
  x<-rle(fixes$ID)
  events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=1))
  events=events[events[,2]>=dur,]
  include<-match(fixes[,1],events[,1])
  J<-which(is.na(include)==0)
  fixes=fixes[J,]
  }
fixes=fixes[fixes$CV>=thresh,]
tmp=table(factor(fixes$Lon2,levels=0:359),factor(fixes$Lat2,levels=-90:89),fixes$Month)
systems[,,y,]=tmp
print(mean(systems[,,y,],na.rm=T))
} # End year loop

## Write a netcdf file with cyclones - so can read later

#dimX<-dim.def.ncdf("lon","degrees_E",lon)
#dimY<-dim.def.ncdf("lat","degrees_N",lat)
#dimT1<-dim.def.ncdf("year","years",years)
#dimT2<-dim.def.ncdf("month","months",1:12)

dimX<-ncdim_def("lon","degrees_E",lon)
dimY<-ncdim_def("lat","degrees_N",lat)
dimT1<-ncdim_def("year","years",years)
dimT2<-ncdim_def("month","months",months)

fillvalue <- 1e32
cyc_def <- ncvar_def("systems","count",list(dimX,dimY,dimT1,dimT2),fillvalue,paste("Number of",type,"pressure systems during each month for each location in the Australian region"),prec="single")

# create netCDF file and put arrays
if(is.na(dur)) ncfname <- paste(dir,"/count_",type,"s","_cv",thresh,".nc",sep="") else ncfname<-paste(dir,"/count_",type,"s","_cv",thresh,"_D",dur,".nc",sep="")
ncout <- nc_create(ncfname,cyc_def) #force_v4=T)

# put variables
ncvar_put(ncout,cyc_def,systems)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"year","axis","T1")
ncatt_put(ncout,"month","axis","T2")

nc_close(ncout)

} # End function

#extract_counts(1950,2016,type="low",dir="gcyc_out",thresh=0.25,dur=3)
#extract_counts(1990,2012,type="low",dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj240_rad5cv0.15/",thresh=0.15)

#extract_counts(1980,2016,type="high",dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad5cv0.075/",thresh=0.15,dur=3)
extract_counts(1950,2016,type="high",dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_highs_rad10cv0.075/",thresh=0.075,dur=3)

