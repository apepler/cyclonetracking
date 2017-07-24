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
ECLs<-array(0,c(length(years),length(months),length(members),4))

lat=seq(-45,-9,2)
lon=seq(105,165,2)
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
fixes=cbind(fixes,matrix(0,dim(fixes)[1],2))

### Part 1 - identify ECLs

  I<-which(fixes[,6]>=149 & fixes[,6]<=161 & fixes[,7]<(-37) & fixes[,7]>=-41)
  fixes[I,11]<-1
  I<-which(fixes[,6]>=(149+(37+fixes[,7])/2) & fixes[,6]<=161 & fixes[,7]<(-31) & fixes[,7]>=-37)
  fixes[I,11]<-1
  I<-which(fixes[,6]>=152 & fixes[,6]<=161 & fixes[,7]<=(-24) & fixes[,7]>=-31)
  fixes[I,11]<-1
  fixes[,12]=floor(fixes[,3]/100)%%100 #Month

  months2=seq(months[m],by=1,length.out=4)
  months2[months2>12]=months2[months2>12]-12

  for(m2 in 1:4)
  {
  I=which(fixes[,11]==1 & fixes[,9]>=thresh & fixes[,12]==months2[m2])
  ECLs[y,m,e,m2]=length(I)
  }

#### Part 2 - add all cyclones to a grid

for(i in 1:length(lon))
  for(j in 1:length(lat))
    for(m2 in 1:length(months2)) 
    {
     I=which(fixes[,6]>=lon[i]-1 & fixes[,6]<lon[i]+1 &
               fixes[,7]>=lat[j]-1 & fixes[,7]<lat[j]+1 &
               fixes[,9]>=thresh & fixes[,12]==months2[m2]) 
     cyclones[i,j,y,m,e,m2]=length(I)
    }

} # End time loop

## Next step - write ECL file

dimT1<-ncdim_def("year","years",years)
dimT2<-ncdim_def("month","months",months)
dimN<-ncdim_def("members","count",1:length(members))
dimM<-ncdim_def("leadtime","months",0:3)

fillvalue <- 1e32
ECL_def <- ncvar_def("ECLs","count",list(dimT1,dimT2,dimN,dimM),fillvalue,"Number of ECLs during each month & member",prec="single")

# create netCDF file and put arrays
ncfname <- paste(dir,"/ECLs_4monthlead_day",inday,"_cv",thresh,".nc",sep="")
ncout <- nc_create(ncfname,ECL_def,force_v4=T)

# put variables
ncvar_put(ncout,ECL_def,ECLs)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"year","axis","T1")
ncatt_put(ncout,"month","axis","T2")
ncatt_put(ncout,"members","axis","N")
ncatt_put(ncout,"leadtime","axis","M")

nc_close(ncout)


## Write a netcdf file with cyclones - so can read later

dimX<-ncdim_def("lon","degrees_E",lon)
dimY<-ncdim_def("lat","degrees_N",lat)
cyc_def <- ncvar_def("cyclones","count",list(dimX,dimY,dimT1,dimT2,dimN,dimM),fillvalue,"Number of cyclones during each month & member for each location in the Australian region",prec="single")

# create netCDF file and put arrays
ncfname <- paste(dir,"/austcyclones_4monthlead_day",inday,"_cv",thresh,".nc",sep="")
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

#This is a wrapper function for colorRampPalette. It allows for the
#definition of the number of intermediate colors between the main colors.
#Using this option one can stretch out colors that should predominate
#the palette spectrum. Additional arguments of colorRampPalette can also
#be added regarding the type and bias of the subsequent interpolation.
color.palette <- function(steps, n.steps.between=NULL, ...){
  
  if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
  if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
  
  fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
  RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
  RGB[,fill.steps] <- col2rgb(steps)
  
  for(i in which(n.steps.between>0)){
    col.start=RGB[,fill.steps[i]]
    col.end=RGB[,fill.steps[i+1]]
    for(j in seq(3)){
      vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
      RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
    }
  }
  
  new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
  pal <- colorRampPalette(new.steps, ...)
  return(pal)
}

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  poly <- vector(mode="list", length(cols))
  for(i in seq(poly)){
    poly[[i]] <- c(i, i+1, i+1, i)
  }
  plot(1,1,t="n",ylim=c(0,1), xlim=c(1,length(cols)), axes=F,xlab='',ylab='')  
  for(i in seq(poly)){
    polygon(poly[[i]], c(0,0,1,1), col=cols[i], border=NA)
  }
  axis(1, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}

ACCESS_ECLs_4month(1990,2012,1,thresh=0.15)
