       program getmonthyear

cc     routine that takes a month (or season) and year input and outputs the appropriate
cc     character strings for start and end (month/season and year) periods

cc     language: fortran90
cc     author: K. Braganza

cc     last modified: Nov 2006


cc **************************************************************************************************
cc SUBROUTINE TO CREATE AUSMAP OUTPUT

       implicit none

       integer
     $ xmonth,imonth1,imonth2
     $,xyear,iyear1,iyear2
     $,cday1,cday2
     $,season_flag

       character
     $ cdate1*8
     $,cdate2*8

       character
     $ cyear*8
     $,cmonth*4
     $,cseason*4

       character*10 
     $ amonth(12)

       integer
     $ lenmonth(12)

       double precision
     $ test1,test2,test3

       logical
     $ divisible_by_4
     $,divisible_by_100
     $,divisible_by_400
     $,days_30
     $,leap_year
     $,february
     $,djf,mam,jja,son,nws,sws,ann
     $,season

cc **************************************************************************************************
cc     Set the timestamp header from ausmap conventions

       amonth = (/ 'January   ','February  ','March     ',
     $           'April     ','May       ','June      ',
     $           'July      ','August    ','September ',
     $           'October   ','November  ','December  '/)

       lenmonth = (/ 7,8,5,5,3,4,4,6,9,7,8,8 /)

       call getarg(1,cyear)
       call getarg(2,cmonth)
       call getarg(3,cseason)

       read(cyear,*,err=50)xyear 
       read(cmonth,*,err=50)xmonth 
       read(cseason,*,err=50)season_flag 

       season = .FALSE.

c      set the defaults

       days_30 = .FALSE.
       cday1 = 1
       cday2 = 31

       leap_year = .FALSE.

       iyear1=xyear
       iyear2=xyear

       imonth1=xmonth
       imonth2=xmonth

c      set seasonal file conventions

       season = season_flag.eq.1

       if(season) then   ! if season flag is on use seasonal dates

         djf = xmonth.eq.1
         mam = xmonth.eq.2
         jja = xmonth.eq.3
         son = xmonth.eq.4
         nws = xmonth.eq.5         ! Northern Wet Season
         sws = xmonth.eq.6         ! Southern Wet Season
         ann = xmonth.eq.7         ! Annual

          if(djf)then
           iyear1=xyear-1
           iyear2=xyear
           imonth1=12
           imonth2=2
          elseif(mam)then
           imonth1=3
           imonth2=5
          elseif(jja)then
           imonth1=6
           imonth2=8
          elseif(son)then
           imonth1=9
           imonth2=11
          elseif(nws)then
           iyear1=xyear-1
           iyear2=xyear
           imonth1=10
           imonth2=4
          elseif(sws)then
           imonth1=4
           imonth2=11
          elseif(ann)then
           imonth1=1
           imonth2=12
          endif

       endif   ! seasonal condition


c      set the days of the month 30.or.31

       days_30 = 
     $ (imonth2.eq.4).or.(imonth2.eq.6)
     $  .or.(imonth2.eq.9).or.(imonth2.eq.11)

       if(days_30)then
          cday2 = 30
       endif

c      set the days of the month for february leap.or.non_leap

       !open(unit=21,file='leap_years.out',form='formatted')

       test1=float(iyear2)/4.0
       test2=float(iyear2)/100.0
       test3=float(iyear2)/400.0

       divisible_by_4 = test1.eq.NINT(test1)
       divisible_by_100 = test2.eq.NINT(test2)
       divisible_by_400 = test3.eq.NINT(test3)
    
       if(divisible_by_4)then
        leap_year = .TRUE.
        !write(21,*)iyear2,test1,NINT(test1),'LEAP YEAR!'
         if(divisible_by_100)then
          leap_year = .FALSE.
          !write(21,*)iyear2,test2,NINT(test2),'not a leap year'
         endif
         if(divisible_by_400)then
          leap_year = .TRUE.
          !write(21,*)iyear2,test3,NINT(test3),'LEAP YEAR'
         endif
       else
        leap_year = .FALSE.
        !write(21,*)iyear2,'not a leap year'
       endif
        
       february = imonth2.eq.2

       if(february)then
         if(leap_year)then
           cday2 = 29
         else
           cday2 = 28
         endif
       endif

       write(cdate1,'(i4.4,i2.2,i2.2)')iyear1,imonth1,cday1
       write(cdate2,'(i4.4,i2.2,i2.2)')iyear2,imonth2,cday2

       
       open(unit=11,file='getmonthyear.out')

       write(11,1100)'STARTDATE ',cdate1
       write(11,1100)'ENDDATE   ',cdate2
        
1100   format(A10,A8)

       if(.not.season)then
       write(11,*) 
     $,amonth(imonth1)(1:lenmonth(imonth1))
     $,iyear1
       endif

       if(season)then
        if(iyear1.eq.iyear2)then
       write(11,*) 
     $ cday1
     $,amonth(imonth1)(1:lenmonth(imonth1))
     $,' to'
     $,cday2
     $,amonth(imonth2)(1:lenmonth(imonth2))
     $,iyear1
        else
       write(11,*) 
     $ cday1
     $,amonth(imonth1)(1:lenmonth(imonth1))
     $,iyear1            
     $,'to'
     $,cday2
     $,amonth(imonth2)(1:lenmonth(imonth2))
     $,iyear2
        endif
       endif

       close(unit=11)

       write(6,*) 'GETMONTHYEAR NORMAL TERMINATION' 
      
       stop

50     continue

       write(6,*) 'GETMONTHYEAR INPUT ERROR' 
      
       stop

       end
