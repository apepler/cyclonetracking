## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

library(ncdf4)

ACCESS_ECLs_4month<-function(year1,year2,inday,dir="../gcyc_out/access-s1/proj100_rad5cv0.15/",thresh=0.15)
{
years=seq(year1,year2,1)
#months=1:12
months=5 ## For the test version
members=paste("e",sprintf(0:10,fmt="%2.2d"),sep="")

lat=seq(-89.5,89.5,1)
lon=seq(0.5,359.5,1)

members=paste("e",sprintf(1:11,fmt="%2.2d"),sep="")

cyclones<-array(0,c(length(lon),length(lat),length(years),length(months),length(members),4))

for(y in 1:length(years))
for(m in 1:length(months))
for(e in 1:length(members))
{
indate=paste(years[y],sprintf("%02d",months[m]),sprintf("%02d",inday),sep="")
fname=paste(dir,members[e],"/tracks_",indate,"_4month.dat",sep="")
print(fname)
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=floor(fixes$Lat)
fixes$Lon2=floor(fixes$Lon)%%360

fixes=fixes[fixes$CV>=thresh,]
tmp=table(factor(fixes$Lon2,levels=0:359),factor(fixes$Lat2,levels=-90:89),fixes$Month)
cyclones[,,y,m,e,]=tmp


} # End time loop

## Next step - write ECL file

dimT1<-ncdim_def("year","years",years)
dimT2<-ncdim_def("month","months",months)
dimN<-ncdim_def("members","count",1:length(members))
dimM<-ncdim_def("leadtime","months",0:3)

fillvalue <- 1e32

## Write a netcdf file with cyclones - so can read later

dimX<-ncdim_def("lon","degrees_E",lon)
dimY<-ncdim_def("lat","degrees_N",lat)
cyc_def <- ncvar_def("cyclones","count",list(dimX,dimY,dimT1,dimT2,dimN,dimM),fillvalue,"Number of cyclones during each month & member for each location",prec="single")

# create netCDF file and put arrays
ncfname <- paste(dir,"/globalcyclones_4monthlead_day",inday,"_cv",thresh,".nc",sep="")
ncout <- nc_create(ncfname,cyc_def,force_v4=T)

# put variables
ncvar_put(ncout,cyc_def,cyclones)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"year","axis","T1")
ncatt_put(ncout,"month","axis","T2")
ncatt_put(ncout,"members","axis","N")
ncatt_put(ncout,"leadtime","axis","M")

nc_close(ncout)


} # End function


ACCESS_ECLs_4month(1990,2012,1,dir="../gcyc_out/access-s1/proj240_rad5cv0.15/",thresh=0.15)
