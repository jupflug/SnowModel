c outputs_user.f

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine OUTPUTS_USER(nx,ny,iter,Tair_grid,rh_grid,
     &  uwind_grid,vwind_grid,windspd_grid,winddir_grid,
     &  Qsi_grid,Qli_grid,prec_grid,Tsfc,Qle,Qh,Qe,Qc,Qm,Qf,
     &  e_balance,snow_depth,xro_snow,swe_depth,ro_nsnow,
     &  runoff,rain,sprec,sum_prec,sum_runoff,w_balance,
     &  snow_d,topo_land,wbal_qsubl,sum_sprec,wbal_salt,
     &  wbal_susp,ro_snow_grid,sum_Qcs,canopy_int,Qcs,
     &  iyear,imonth,iday,xhour,undef,deltax,xmn,ymn,
     &  wbal_subgrid,canopy_unload,sum_qsubl,sum_trans,
     &  sum_unload,sum_glacmelt,glacier_melt,swemelt,
     &  iprint_inc,sfc_pressure,sum_swemelt,albedo,
     &  icorr_factor_loop,swesublim,vegtype,iter_start,
     &  seaice_run,print_inc)

c This subroutine is available to provide user-defined outputs.
c   These might be special-case situations, like just writing out
c   data at the end of every day, writing out a few grid cells,
c   saving each data arrays to individual files, etc.

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,iter,max_poly,num_poly,iyear,imonth,iday,
     &  iprint_inc,icorr_factor_loop,iter_start,icorr_loop_new,
     &  individual_files,n_vars,k,irec

      real Tair_grid(nx_max,ny_max),rh_grid(nx_max,ny_max),
     &  uwind_grid(nx_max,ny_max),vwind_grid(nx_max,ny_max),
     &  windspd_grid(nx_max,ny_max),winddir_grid(nx_max,ny_max),
     &  Qsi_grid(nx_max,ny_max),Qli_grid(nx_max,ny_max),
     &  prec_grid(nx_max,ny_max),Tsfc(nx_max,ny_max),
     &  Qle(nx_max,ny_max),Qh(nx_max,ny_max),Qe(nx_max,ny_max),
     &  Qc(nx_max,ny_max),Qm(nx_max,ny_max),Qf(nx_max,ny_max),
     &  e_balance(nx_max,ny_max),snow_depth(nx_max,ny_max),
     &  xro_snow(nx_max,ny_max),swe_depth(nx_max,ny_max),
     &  ro_nsnow(nx_max,ny_max),runoff(nx_max,ny_max),
     &  rain(nx_max,ny_max),sprec(nx_max,ny_max),
     &  sum_prec(nx_max,ny_max),sum_runoff(nx_max,ny_max),
     &  w_balance(nx_max,ny_max),snow_d(nx_max,ny_max),
     &  topo_land(nx_max,ny_max),wbal_qsubl(nx_max,ny_max),
     &  sum_sprec(nx_max,ny_max),wbal_salt(nx_max,ny_max),
     &  wbal_susp(nx_max,ny_max),ro_snow_grid(nx_max,ny_max),
     &  sum_Qcs(nx_max,ny_max),canopy_int(nx_max,ny_max),
     &  Qcs(nx_max,ny_max),wbal_subgrid(nx_max,ny_max),
     &  canopy_unload(nx_max,ny_max),sum_qsubl(nx_max,ny_max),
     &  sum_trans(nx_max,ny_max),glacier_melt(nx_max,ny_max),
     &  sum_unload(nx_max,ny_max),sum_glacmelt(nx_max,ny_max),
     &  swemelt(nx_max,ny_max),sfc_pressure(nx_max,ny_max),
     &  sum_swemelt(nx_max,ny_max),swesublim(nx_max,ny_max),
     &  vegtype(nx_max,ny_max),albedo(nx_max,ny_max)

      real undef,xhour,deltax,pi,rad2deg,seaice_run,print_inc
      double precision xmn,ymn

      real uwnd(nx_max,ny_max)
      real vwnd(nx_max,ny_max)

c Define the output variable data block.
      parameter (n_vars=20)
      real vars(nx_max,ny_max,n_vars)

      character*1 c_var(n_vars)
      character*4 c_outvars(n_vars)

      data c_outvars /'tair','relh','wspd','qsin','qlin',
     &                'qlem','albd','wdir','prec','rpre',
     &                'spre','smlt','ssub','roff','glmt',
     &                'snod','sden','swed','sspr','ssmt'/
      
      character path1*(*) 
      character path2*(*) 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c BEGIN USER EDIT SECTION.

c Define which variables you want to save, whether they will be
c   output at every time step or some time-step increment, which
c   directories the data will be put in, etc.

c Define the output file locations (paths).
      parameter (path1 =
     &  'outputs/wo_assim/') 
      parameter (path2 =
     &  'outputs/wi_assim/') 

c Write a seperate file for each variable (individual_files = 1).
c   No other option has been implemented here.
      individual_files = 1

c Define the number of time steps you are going to sum or average
c   over.  If you want to output data at every model time step, set
c   print_inc = 1.0.  For run with an hourly time step and data
c   writes once a day, print_inc = 24.0.  For a run with 3-hourly
c   time steps and data writes once a day, print_inc = 8.0.
      print_inc = 24.0
c     print_inc = 1.0
c     print_inc = 8.0

c Define the variables you want to save.  The following are the
c   variables this subroutine is currently set up to output.
c   Listed are the output variable name and the corresponding model
c   variable name.

c VALUES AVERAGED OVER THE PERIOD.
c    1   tair(i,j) = Tair_grid(i,j) - 273.16
c    2   relh(i,j) = rh_grid(i,j)
c    3   wspd(i,j) = windspd_grid(i,j)
c    4   qsin(i,j) = Qsi_grid(i,j)
c    5   qlin(i,j) = Qli_grid(i,j)
c    6   qlem(i,j) = Qle(i,j)
c    7   albd(i,j) = albedo(i,j)
c    8   wdir(i,j) = from uwind_grid(i,j) and vwind_grid(i,j)

c VALUES SUMMED OVER THE PERIOD.
c    9   prec(i,j) = prec_grid(i,j)
c   10   rpre(i,j) = rain(i,j)
c   11   spre(i,j) = sprec(i,j)
c   12   smlt(i,j) = swemelt(i,j)
c   13   ssub(i,j) = swesublim(i,j)
c   14   roff(i,j) = runoff(i,j)
c   15   glmt(i,j) = glacier_melt(i,j)

c VALUES SAVED AT THE END OF THE PERIOD.
c   16   snod(i,j) = snow_depth(i,j)
c   17   sden(i,j) = xro_snow(i,j)
c   18   swed(i,j) = swe_depth(i,j)
c   19   sspr(i,j) = sum_sprec(i,j)
c   20   ssmt(i,j) = sum_swemelt(i,j)

c Define which variables you want to save by placing a yes = 'y' or
c   no = 'n' in front of the variable number.

c VALUES AVERAGED OVER THE PERIOD.
c VALUES AVERAGED OVER THE PERIOD.
c 1 = tair(i,j) = Tair_grid(i,j) - 273.16
      c_var(1)  = 'n'

c 2 = relh(i,j) = rh_grid(i,j)
      c_var(2)  = 'n'

c 3 = wspd(i,j) = windspd_grid(i,j)
      c_var(3)  = 'n'

c 4 = qsin(i,j) = Qsi_grid(i,j)
      c_var(4)  = 'n'

c 5 = qlin(i,j) = Qli_grid(i,j)
      c_var(5)  = 'n'

c 6 = qlem(i,j) = Qle(i,j)
      c_var(6)  = 'n'

c 7 = albd(i,j) = albedo(i,j)
      c_var(7)  = 'n'

c 8 = wdir(i,j) = from uwind_grid(i,j) and vwind_grid(i,j)
      c_var(8)  = 'n'

c VALUES SUMMED OVER THE PERIOD.
c VALUES SUMMED OVER THE PERIOD.
c  9 = prec(i,j) = prec_grid(i,j)
      c_var(9)  = 'n'

c 10 = rpre(i,j) = rain(i,j)
      c_var(10) = 'n'

c 11 = spre(i,j) = sprec(i,j)
      c_var(11) = 'n'

c 12 = smlt(i,j) = swemelt(i,j)
      c_var(12) = 'n'

c 13 = ssub(i,j) = swesublim(i,j)
      c_var(13) = 'n'

c 14 = roff(i,j) = runoff(i,j)
      c_var(14) = 'n'

c 15 = glmt(i,j) = glacier_melt(i,j)
      c_var(15) = 'n'

c VALUES SAVED AT THE END OF THE PERIOD.
c VALUES SAVED AT THE END OF THE PERIOD.
c 16 = snod(i,j) = snow_depth(i,j)
      c_var(16) = 'y'

c 17 = sden(i,j) = xro_snow(i,j)
      c_var(17) = 'n'

c 18 = swed(i,j) = swe_depth(i,j)
      c_var(18) = 'y'

c 19 = sspr(i,j) = sum_sprec(i,j)
      c_var(19) = 'n'

c 20 = ssmt(i,j) = sum_swemelt(i,j)
      c_var(20) = 'n'

c Note that this data output implementation is currently configured
c   to mask out the ocean points (vegtype.eq.24.0) if this is a
c   land run (seaice_run = 0.0); mask out all land points (vegtype.
c   ne.24.0) if this is an ocean/sea ice run (seaice_run = 1.0); and
c   to not mask out anything if this is a combined land and sea ice
c   run (seaice_run = 2.0).

c END USER EDIT SECTION.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Define the constants used in the wind-direction averaging.
      pi = 2.0 * acos(0.0)
      rad2deg = 180.0 / pi

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Use individual output files for each variable.
      if (individual_files.eq.1) then

c Open individual output files for each variable.
        if (iter.eq.iter_start) then
          if (icorr_factor_loop.eq.1) then
            do k=1,n_vars
              if (c_var(k).eq.'y') then
                open (220+k,file=path1//c_outvars(k)//'.gdat',
     &            form='unformatted',access='direct',recl=4*nx*ny)
              endif
            enddo
          endif

          if (icorr_factor_loop.eq.2) then
            do k=1,n_vars
              if (c_var(k).eq.'y') then
                open (320+k,file=path2//c_outvars(k)//'.gdat',
     &            form='unformatted',access='direct',recl=4*nx*ny)
              endif
            enddo
          endif
        endif

        if (iter.eq.iter_start) then
c Initialize the averaging and summing arrays.
          do j=1,ny
            do i=1,nx
              do k=1,n_vars
                vars(i,j,k) = 0.0
              enddo
              uwnd(i,j) = 0.0
              vwnd(i,j) = 0.0
            enddo
          enddo
        endif

c Perform the avaraging, summing, etc.
        do j=1,ny
          do i=1,nx
c Values averaged over the period.
            vars(i,j,1) = vars(i,j,1) + (Tair_grid(i,j) - 273.16) /
     &        print_inc
            vars(i,j,2) = vars(i,j,2) + rh_grid(i,j) / print_inc
            vars(i,j,3) = vars(i,j,3) + windspd_grid(i,j) / print_inc
            vars(i,j,4) = vars(i,j,4) + Qsi_grid(i,j) / print_inc
            vars(i,j,5) = vars(i,j,5) + Qli_grid(i,j) / print_inc
            vars(i,j,6) = vars(i,j,6) + Qle(i,j) / print_inc
            vars(i,j,7) = vars(i,j,7) + albedo(i,j) / print_inc

            uwnd(i,j) = uwnd(i,j) + uwind_grid(i,j) / print_inc
            vwnd(i,j) = vwnd(i,j) + vwind_grid(i,j) / print_inc

c Some compilers do not allow both u and v to be 0.0 in
c   the atan2 computation.
            if (abs(uwnd(i,j)).lt.1e-10) uwnd(i,j) = 1e-10

            vars(i,j,8) = rad2deg * atan2(uwnd(i,j),vwnd(i,j))
            if (vars(i,j,8).ge.180.0) then
              vars(i,j,8) = vars(i,j,8) - 180.0
            else
              vars(i,j,8) = vars(i,j,8) + 180.0
            endif

c Values summed over the period.
            vars(i,j,9) = vars(i,j,9) + prec_grid(i,j)
            vars(i,j,10) = vars(i,j,10) + rain(i,j)
            vars(i,j,11) = vars(i,j,11) + sprec(i,j)
            vars(i,j,12) = vars(i,j,12) + swemelt(i,j)
            vars(i,j,13) = vars(i,j,13) + swesublim(i,j)
            vars(i,j,14) = vars(i,j,14) + runoff(i,j)
            vars(i,j,15) = vars(i,j,15) + glacier_melt(i,j)
c Values saved at the end of the day.
            vars(i,j,16) = snow_depth(i,j)
            vars(i,j,17) = xro_snow(i,j)
            vars(i,j,18) = swe_depth(i,j)
            vars(i,j,19) = sum_sprec(i,j)
            vars(i,j,20) = sum_swemelt(i,j)
          enddo
        enddo

c Check to see whether this is the data-write time step.
        if (mod(iter,nint(print_inc)).eq.0) then

c Mask out the ocean points (vegtype.eq.24.0) if this is
c   a land run (seaice_run = 0.0).  Mask out all land points
c   (vegtype.ne.24.0) if this is an ocean/sea ice run
c   (seaice_run = 1.0).  Do not mask out anything if this is
c   a combined land and sea ice run (seaice_run = 2.0).
          if (seaice_run.eq.0.0) then
            do j=1,ny
              do i=1,nx
                if (vegtype(i,j).eq.24.0) then
                  do k=1,n_vars
                    vars(i,j,k) = undef
                  enddo
                endif
              enddo
            enddo
          elseif (seaice_run.eq.1.0) then
            do j=1,ny
              do i=1,nx
                if (vegtype(i,j).ne.24.0) then
                  do k=1,n_vars
                    vars(i,j,k) = undef
                  enddo
                endif
              enddo
            enddo
          endif

c Write out the data.
          irec = iter / nint(print_inc)
          if (icorr_factor_loop.eq.1) then
            do k=1,n_vars
              if (c_var(k).eq.'y') then
                write (220+k,rec=irec) ((vars(i,j,k),i=1,nx),j=1,ny)
              endif
            enddo
          elseif (icorr_factor_loop.eq.2) then
            do k=1,n_vars
              if (c_var(k).eq.'y') then
                write (320+k,rec=irec) ((vars(i,j,k),i=1,nx),j=1,ny)
              endif
            enddo
          endif

c Reinitialize the averaging and summing arrays.
          do j=1,ny
            do i=1,nx
              do k=1,n_vars
                vars(i,j,k) = 0.0
              enddo
              uwnd(i,j) = 0.0
              vwnd(i,j) = 0.0
            enddo
          enddo
        endif

c Use more than one variable in an output file.
      else

        print *,'Use more than one variable in an output file:'
        print *,'  THIS HAS NOT BEEN IMPLEMENTED YET'

      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c THIS IS AN EXAMPLE OF SAVING DATA IN ASCII/TEXT FORMAT.

c Save the swe_depth data at the end of every day.
c     if (mod(iter,iprint_inc).eq.0) then
c       call ascii_outputs_1(nx,ny,iyear,imonth,iday,xhour,undef,
c    &    swe_depth,deltax,xmn,ymn)
c     endif

c Save the Tair data every hour for the first day of the simulation.
c     if (iter.le.24) then
c       call ascii_outputs_2(nx,ny,iyear,imonth,iday,xhour,undef,
c    &    Tair_grid,deltax,xmn,ymn)
c     endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ascii_outputs_1(nx,ny,iyear,imonth,iday,xhour,undef,
     &  swe_depth,deltax,xmn,ymn)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,iyear,imonth,iday,ihour
      real undef,xhour,deltax
      double precision xmn,ymn
      real swe_depth(nx_max,ny_max)

      character*18 name1
      character*1 dot
      character*4 name2
      character*4 yyyy
      character*2 mm
      character*2 dd
      character*2 hh
      character*35 fname
      character*40 form

      name1 = 'outputs/swe_depth_'
      name2 = '.asc'
      dot = '.'

      ihour = nint(xhour)
      write(yyyy,'(i4.4)') iyear
      write(mm,'(i2.2)') imonth
      write(dd,'(i2.2)') iday
      write(hh,'(i2.2)') ihour
      fname = name1//yyyy//dot//mm//dot//dd//dot//hh//name2

      open (23,file=fname)

      write (23,*) 'ncols        ',nx
      write (23,*) 'nrows        ',ny
      write (23,*) 'xllcorner    ',xmn
      write (23,*) 'yllcorner    ',ymn
      write (23,*) 'cellsize     ',deltax
      write (23,*) 'NODATA_value ',undef

c Define the output format.  The following will produce nx columns
c   in the ascii output array (like assumed in an ARC/INFO GRID
c   ascii data file.  nx is getting written to the i5 space.  The
c   output format for the data values will be f12.4 (this is what
c   you might want to modify, depending on the data units/values).
      write (form,90) nx
  90  format ('(',i5,'f12.4)')

      do j=ny,1,-1
        write (23,form) (swe_depth(i,j),i=1,nx)
      enddo

      close (23)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ascii_outputs_2(nx,ny,iyear,imonth,iday,xhour,undef,
     &  Tair_grid,deltax,xmn,ymn)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,iyear,imonth,iday,ihour
      real undef,xhour,deltax
      double precision xmn,ymn
      real Tair_grid(nx_max,ny_max)

      character*13 name1
      character*1 dot
      character*4 name2
      character*4 yyyy
      character*2 mm
      character*2 dd
      character*2 hh
      character*30 fname
      character*40 form

      name1 = 'outputs/Tair_'
      name2 = '.asc'
      dot = '.'

      ihour = nint(xhour)
      write(yyyy,'(i4.4)') iyear
      write(mm,'(i2.2)') imonth
      write(dd,'(i2.2)') iday
      write(hh,'(i2.2)') ihour
      fname = name1//yyyy//dot//mm//dot//dd//dot//hh//name2

      open (23,file=fname)

      write (23,*) 'ncols        ',nx
      write (23,*) 'nrows        ',ny
      write (23,*) 'xllcorner    ',xmn
      write (23,*) 'yllcorner    ',ymn
      write (23,*) 'cellsize     ',deltax
      write (23,*) 'NODATA_value ',undef

c Define the output format.  The following will produce nx columns
c   in the ascii output array (like assumed in an ARC/INFO GRID
c   ascii data file.  nx is getting written to the i5 space.  The
c   output format for the data values will be f12.4 (this is what
c   you might want to modify, depending on the data units/values).
      write (form,90) nx
  90  format ('(',i5,'f12.4)')

      do j=ny,1,-1
        write (23,form) (Tair_grid(i,j)-273.16,i=1,nx)
      enddo

      close (23)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c THE CODE BELOW WAS USED TO SAVE AVERAGES OVER POLYGONS.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c These parameters and variables control polygon outputs.
c     parameter (max_poly=1002001)
c     real poly(nx_max,ny_max)
c     real tair_poly(max_poly)
c     real rh_poly(max_poly)
c     real wspd_poly(max_poly)
c     real prec_poly(max_poly)
c     real Qsi_poly(max_poly)
c     real Qli_poly(max_poly)
c     real count_poly(max_poly)
c     real elke
c     character*7 ipoly_num
c     character*32 outfname

c     elke = 0.0
c     if (elke.eq.1.0) then
c Perform the preprocessing for the polygon outputs to be used by
c   FASST and SNTHERM.
c       if (iter.eq.1) then
c         open (47,file='polys/fraser.polys.200m.gdat',
c    &      form='unformatted',access='direct',recl=4*nx*ny)

c         read (47,rec=1) ((poly(i,j),i=1,nx),j=1,ny)

c Generate a table describing the number of grid cells in each
c   polygon.
c         do num_poly=1,max_poly
c           count_poly(num_poly) = 0.0
c         enddo
c         do j=1,ny
c           do i=1,nx
c             num_poly = nint(poly(i,j))
c             if (num_poly.ne.undef)
c    &        count_poly(num_poly) = count_poly(num_poly) + 1.0
c           enddo
c         enddo

c Open the polygon output files.
c         do num_poly=1,max_poly
c           if (count_poly(num_poly).gt.0.0) then
c             write(ipoly_num,'(i7.7)') num_poly
c             outfname = 'polys/fraser.polygon-'//ipoly_num//'.dat'
c             open (num_poly+100,file=outfname,form='formatted')
c           endif
c         enddo
c       endif

c Initialize the polygon averaging array.
c       do num_poly=1,max_poly
c         tair_poly(num_poly) = 0.0
c         rh_poly(num_poly) = 0.0
c         wspd_poly(num_poly) = 0.0
c         prec_poly(num_poly) = 0.0
c         Qsi_poly(num_poly) = 0.0
c         Qli_poly(num_poly) = 0.0
c       enddo

c  Calculate the polygon averages be used by FASST and SNTHERM.
c       call ave_poly_var(nx,ny,poly,count_poly,tair_poly,Tair_grid,
c    &    undef,max_poly)

c       call ave_poly_var(nx,ny,poly,count_poly,rh_poly,rh_grid,
c    &    undef,max_poly)

c       call ave_poly_var(nx,ny,poly,count_poly,wspd_poly,windspd_grid,
c    &    undef,max_poly)

c       call ave_poly_var(nx,ny,poly,count_poly,prec_poly,prec_grid,
c    &    undef,max_poly)

c       call ave_poly_var(nx,ny,poly,count_poly,Qsi_poly,Qsi_grid,
c    &    undef,max_poly)

c       call ave_poly_var(nx,ny,poly,count_poly,Qli_poly,Qli_grid,
c    &    undef,max_poly)

c Write the data to individual polygon files.
c       do num_poly=1,max_poly
c         if (count_poly(num_poly).gt.0.0)
c    &      write (num_poly+100,88) iyear,imonth,iday,xhour,
c    &        num_poly,tair_poly(num_poly)-273.16,rh_poly(num_poly),
c    &        wspd_poly(num_poly),prec_poly(num_poly),
c    &        Qsi_poly(num_poly),Qli_poly(num_poly)
c       enddo

c 88    format (i5,i3,i3,f6.2,i8,3f10.3,f12.6,2f10.3)
c     endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     subroutine ave_poly_var(nx,ny,poly,count_poly,ave_poly,var,
c    &  undef,max_poly)

c     implicit none

c     include 'snowmodel.inc'

c     integer i,j,nx,ny,num_poly,max_poly
c     real undef
c     real poly(nx_max,ny_max)
c     real var(nx_max,ny_max)
c     real ave_poly(max_poly)
c     real count_poly(max_poly)

c Calculate the average met forcing value for each polygon.
c     do j=1,ny
c       do i=1,nx
c         num_poly = nint(poly(i,j))
c         if (num_poly.ne.undef)
c    &      ave_poly(num_poly) = ave_poly(num_poly) + var(i,j)
c       enddo
c     enddo

c     do num_poly=1,max_poly
c       if (count_poly(num_poly).gt.0.0) then
c         ave_poly(num_poly) = ave_poly(num_poly)/count_poly(num_poly)
c       endif
c     enddo

c     return
c     end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

