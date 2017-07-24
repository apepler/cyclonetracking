      program gausmake

c----------------------------------------------------------------------
c
c     PROGRAM: GAUSMAKE
c
c     FILE: gausmake.f
c
c     LANGUAGE: f90
c
c     A program that makes the topography file 
c     Takes in an instruction file, but the user will need to correct
c     If variable names aren't lat/lon/hgt
c
c     USAGE:
c
c     gausmake
c
c     EXPLANATION:
c
c
c     The following information is supplied by the instruction file.
c     Which is provided on the command line
c
c     rtype         Reanalysis type (source)
c     gtype         Grid resolution
c     ilon,jlat     Length of longitude/latitude
c     geopt/unpack  Whether the source is in geopotential form 
c                   and needs unpacking (T/F)
c     ofile/infile  Names of output binary & input netcdf files
c
c     HISTORY:
c
c     Created: A.S. Pepler
c     Based on gauscheck (D.A. Jones)
c     LAST CHANGE: 24/07/2017 
c
c----------------------------------------------------------------------

      implicit none
      character*200 infile,instfile
      character*20 ofile
      character*100 header
      character upack,geopt
      character*5 rtype
      character*8 gtype
      integer inet1,ncopn,ncvid,utopen,icode,latv,lonv,ivar
      integer i,j,k,ilon,jlat,alon,ind,ind2
      parameter(alon=600)
      real latitude(alon)
      real latitude2(alon)
      real longitude(alon)
      real longitude2(alon)
      real x(alon,alon)
      real x2(alon,alon)
      real xscale,xoff
      integer*2 work(alon,alon)
      integer*2 miss
      integer start(3)
      integer count(3)

c Open the instruction file
      call get_command_argument(1,instfile)
      open(3,file=instfile,status='old')

      read(3,10)rtype,gtype,ilon,jlat,geopt,upack,ofile,infile
      close(3)
10    format (//a5/a10/i3/i3/a1/a1/a20/a200)

      open(2,file=ofile,form='unformatted') 
      inet1=ncopn(infile,0,icode)

      latv=ncvid(inet1,'lat',icode)
      call ncvgt(inet1,latv,1,jlat,latitude,icode)
      lonv=ncvid(inet1,'lon',icode)
      call ncvgt(inet1,lonv,1,ilon,longitude,icode)  

      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=ilon
      count(2)=jlat
      count(3)=1
      ivar=ncvid(inet1,'hgt',icode)
       
      if (upack.eq.'T') then
        call ncvgt(inet1,ivar,start,count,work,icode)
        call ncagt(inet1,ivar,'scale_factor',xscale,icode)
        call ncagt(inet1,ivar,'add_offset',xoff,icode)
        call ncagt(inet1,ivar,'missing_value',miss,icode)      
        call unpack(work,x,xscale,xoff,miss,ilon,jlat)
      else
        call ncvgt(inet1,ivar,start,count,x,icode)
      endif

      if (geopt.eq.'T') then
        x=x/9.80665
      endif
CC Create the header line      

      write(header,'(A8,22x,A5,5x,I4.4,I2.2,I2.2,1x,I4.4,4x,A8,5x,A10)')
     $"Z",rtype,1980,1,1,0,"m",gtype

c Fix up any issues with longitudes and latitudes

      if (latitude(jlat).lt.latitude(1)) then
        latitude2=latitude(jlat:1:-1)
        x2=x(:,jlat:1:-1)
        latitude=latitude2
        x=x2
      endif

      if (minval(longitude).lt.0) then
c Find where the negative longitudes end
        ind=minloc(longitude, 1, mask=longitude.gt.0)
        ind2=ilon-ind
        longitude2(1:ind2)=longitude((ind+1):ilon)
        longitude2((ind2+1):ilon)=longitude(1:ind)+360
        x2(1:ind2,:)=x((ind+1):ilon,:)
        x2((ind2+1):ilon,:)=x(1:ind,:)
        longitude=longitude2
        x=x2
      endif

CC Write to the file that's currently open    
      write(2)jlat
      write(2)(latitude(j),j=1,jlat)
      write(2)ilon
      write(2)(longitude(i),i=1,ilon)
      write(2)header
      write(2)((x(i,j),i=1,ilon),j=1,jlat)
      
      return
      end


********************************************************************
CC    These two subroutines (unpack and xhour) are borrowed directly from NCEP
CC    http://www.esrl.noaa.gov/psd/data/gridded/readgeneral.F

      subroutine unpack(x,y,xscale,xadd,miss,ilon,jlat)
      parameter (zmiss=1e30)
      real xscale,xadd
      integer*2 x(ilon,jlat),miss
      real y(ilon,jlat)


      do i=1,ilon
         do j=1,jlat
            if(x(i,j).eq.miss)then
               y(i,j)=zmiss
            else
              y(i,j)=(x(i,j)*xscale)+xadd
           endif
         enddo
      enddo

      return
      end



