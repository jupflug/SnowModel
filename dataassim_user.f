c dataassim_user.f

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine DATAASSIM_USER(nx,ny,icorr_factor_index,
     &  corr_factor,max_iter,deltax,deltay,xmn,ymn,nobs_dates,
     &  print_inc,iday_init,imonth_init,iyear_init,depth_assim,
     &  ihrestart_flag)

c Perform the required correction (precipitation and melt) factor
c   calculations.  To run the data assimilation routines requires
c   additional model inputs to be added here (inputs in addition
c   to those in the snowmodel.par file).  Then you must recompile
c   the code before running it.

c This program works when you have data at one or many points
c   (many individual grid cells), for one or more times.  And for
c   data (e.g., average swe) over areas of grid cells; there can be
c   many of these areas, at many different times.

c The input (swe) data file is assumed to be in /snowmodel/data/.
c   See below for the details of how this file is to be configured.

c The data assimilation produces a couple of extra files the user
c   can look at to study what the model did with the input data
c   file (fname_sweobs) that was provided.  First, there is a text
c   file written to /snowmodel/data/ that provides a summary of the
c   calculations that were done (see unit=77 in the code below for
c   what is written to this file).  Second, there is a copy of the
c   precipitation correction surfaces that were generated.  This
c   is also placed in /snowmodel/data/, the file is a GrADS file
c   (called corr_factor.gdat), and there is also a corr_factor.ctl
c   there that can be modified to fit the current simulation.  The
c   data layers in corr_factor.gdat are structured as follows:
c   The number of layers equals the total number of observation
c   times that were assimilated, plus the number of  years in the
c   assimilation.  The order is: the correction surface (cf) for
c   the time between the start of the simulation and observation
c   time 1 in year 1, the cf for the time between obs time 1 in
c   year 1 and obs time 2 in year 1, etc., then a cf==1.0 for the
c   time between the last obs time and the end of year 1 (or the
c   end of the simulation for a 1-year run).  Then this order
c   repeats for the rest of the years of the simulation.  In the
c   GrADS control file (corr_factor.ctl) these layers are assumed
c   to correspond to different times in the data file (although
c   the actual time increment defined in the .ctl file is not
c   really relevant: for example, t=1 corresponds to the cf for
c   obs 1, t=2 is for obs 2, t=3 is for the time between the last
c   obs time and the end of year 1, t=4 is for the obs 1 in year
c   2, etc.).

      implicit none

      include 'snowmodel.inc'

      real deltax,deltay,beta,areas_flag,print_inc
      real sprec_ratio(max_obs_dates),smelt_ratio(max_obs_dates)
      real corr_factor(nx_max,ny_max,max_obs_dates)
      real areas_mask(nx_max,ny_max)

      double precision xmn,ymn
      double precision xstn(nx_max*ny_max),ystn(nx_max*ny_max)

      integer icorr_factor_index(max_time_steps)
      integer iobs_rec(max_obs_dates)
      integer nobs_dates,nx,ny,max_iter,local_assim_flag,iday_init,
     &  imonth_init,iyear_init,nyear,nyears,nobs_total,nobs_total_cfi,
     &  ihrestart_flag
      integer depth_assim

      character*80 fname_swed,fname_sspr,fname_ssmt,fname_sden
      character*80 fname_sweobs
      character*80 fname_sweobs_barnes_mask

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c BEGIN USER EDIT SECTION.
c BEGIN USER EDIT SECTION.

c Define how many years there are in this simulation.
      nyears = 1

c A single data file describes the observation information that
c   will be used in the data assimilation.  This file contains the
c   following information in the following format.  The id can be
c   any number, x and y are easting and northing in m, and swe is
c   in m).

c   total_number_of_observation_dates_for_this_year
c   iyr imo idy (for this observation date)
c   number_of_stations_for_this_observation_date
c   id x y swe
c   id x y swe
c   id x y swe
c   iyr imo idy (for this observation date)
c   number_of_stations_for_this_observation_date
c   id x y swe
c   id x y swe

c For example:

c   2
c   2014 3 15
c   3
c   101 3456.7 23677.4 0.42
c   102 3556.3 25079.3 0.52
c   103 3106.2 29089.3 0.59
c   2014 4 1
c   2
c   101 3456.7 23677.4 0.48
c   103 3106.2 29089.3 0.62

c Then this repeats for each year of the assimilation (the input
c   file looks like the above for a single-year run, and for
c   multi-year runs the data for each following year is just
c   stacked on top of the previous year (like the example below
c   for a two-year assimilation run).  The code requires a
c     "total_number_of_observation_dates_for_this_year"
c   line for each year.  If you have a year with no data, this
c   can be set to be 0 (zero).  If this is 0, the code sets the
c   correction factor to equal 1.0 for that simulation year, and
c   no adjustments are made to the precipitation for that year.

c   2
c   2014 3 15
c   3
c   101 3456.7 23677.4 0.42
c   102 3556.3 25079.3 0.52
c   103 3106.2 29089.3 0.59
c   2014 4 1
c   2
c   101 3456.7 23677.4 0.48
c   103 3106.2 29089.3 0.62
c   1
c   2015 3 25
c   2
c   101 3456.7 23677.4 0.23
c   102 3556.3 25079.3 0.32

c Provide the name of the data file that contains the observed swe
c   information (as described above).
c     fname_sweobs = 'data/Full_pkSWE.txt'
      fname_sweobs = 'data/SDVEl_X276562.5_pkSWEbl.txt'

c Define the file names of the swe depth (swed), annual summed snow
c   precipitation (sspr), and annual summed snowmelt (ssmt) outputs
c   from the first iteration of the data assimilation run.  In this
c   implementation of the data assimilation code, I have assumed
c   that the output files are those created by outputs_user.f,
c   where there is an individual file for each variable.
      fname_swed = 'outputs/wo_assim/swed.gdat'
      fname_sspr = 'outputs/wo_assim/sspr.gdat'
      fname_ssmt = 'outputs/wo_assim/ssmt.gdat'
      fname_sden = 'outputs/wo_assim/sden.gdat'

c THE PARAMETERS BELOW ARE RARELY CHANGED, UNLESS YOU ARE DOING AN
c   AREAS ASSIMILATION (INSTEAD OF ASSIMILATING POINT DATA).

c Beta controls the interpolation distance weights.  Beta = 1.0
c   will give you a very smooth field, and correction factor
c   distributions that may not produce swe's that exactly match
c   the observations.  Beta << 1.0 will give you correction factor
c   fields that go right through the data.  If you just have one
c   data point/area, beta is not used.
      beta = 1.0
c     beta = 0.1
c     beta = 0.5

c Define whether this simulation will be processing areas (data
c   within groups of grid cells: areas_flag = 1.0), or points
c   (single grid cells: areas_flag = 0.0).  Note that if you have
c   a combination of areas and points, you have to use the areas
c   option and treat each point like a single-grid-cell (small)
c   area.
      areas_flag = 0.0
c     areas_flag = 1.0

c If this is an areas simulation, open and read in the areas mask
c   data.  Note that here I assume that the area mask is a nx by ny
c   file with undef values everywhere except at the area 'stations'.
c   And that each 'station' area is given a 1.0, 2.0, etc. that
c   corresponds to the order of the station listing in the 'station'
c   data input file (the first 'station' listed has mask value = 1.0,
c   the second listed has mask value = 2.0, etc.
      if (areas_flag.eq.1.0) then
c       open(63,file='sweobs/nea.obsmask.100m.2009.gdat',
c    &    form='unformatted',access='direct',recl=4*nx*ny)
c       read(63,rec=1) ((areas_mask(i,j),i=1,nx),j=1,ny)
c If you have two masks for two different observation dates, then
c   do something like the following.
c       open(63,file='data/zack_obs_mask.gdat',
c    &    form='unformatted',access='direct',recl=4*nx*ny)
c       read(63,rec=1) ((areas_mask(i,j,1),i=1,nx),j=1,ny)
c       read(63,rec=2) ((areas_mask(i,j,2),i=1,nx),j=1,ny)
      endif

c Define whether this simulation is going to restrict the
c   assimilation influence to some local area surrounding each
c   data point that is assimilated.  This was implemented for
c   ANWR simulations where we only had observations in a corner
c   of the simulation domain and I didn't want the corrections
c   to extend too far outside that local observation area.  So,
c   this is an example of what can be done, and it is not written
c   for general application.  If you want to do something similar,
c   the associated subroutine can be edited for your specific
c   simulation of interest.  For yes, local_assim_flag = 1, for
c   no, local_assim_flag = 0.
      local_assim_flag = 0

c Identify the file that contains the local data assimilation
c   mask.  This is only used if local_assim_flag = 1.
      fname_sweobs_barnes_mask =
     &  '../swe_obs/2014/barnes/obs.gridded.gdat'

c END USER EDIT SECTION.
c END USER EDIT SECTION.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c LOOP THROUGH THE YEARS IN THE SIMULATION.
      do nyear=1,nyears
c Run the data assimilation routines.
        call data_assimilation(nx,ny,deltax,deltay,beta,
     &    areas_flag,sprec_ratio,smelt_ratio,corr_factor,
     &    areas_mask,xmn,ymn,xstn,ystn,nobs_dates,iobs_rec,
     &    local_assim_flag,fname_swed,fname_sspr,fname_ssmt,
     &    fname_sweobs,fname_sweobs_barnes_mask,iday_init,
     &    imonth_init,iyear_init,nyear,nobs_total,depth_assim,
     &    fname_sden)

c Build an array indicating the appropriate correction factor to
c   use at any given time during the simulation.  What this does
c   is define an index array that contains the record number that
c   gets used at every model time step during the second model run
c   loop.  This record number corresponds to the record (krec) of
c   the corr_factor(i,j,krec) array that was generated and saved
c   in the subroutine above.
        call corr_factor_index(nobs_dates,icorr_factor_index,
     &    iobs_rec,max_iter,sprec_ratio,smelt_ratio,print_inc,
     &    iday_init,imonth_init,iyear_init,nyear,nobs_total_cfi,
     &    nyears,ihrestart_flag)

      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine data_assimilation(nx,ny,deltax,deltay,beta,
     &  areas_flag,sprec_ratio,smelt_ratio,corr_factor,
     &  areas_mask,xmn,ymn,xstn,ystn,nobs_dates,iobs_rec,
     &  local_assim_flag,fname_swed,fname_sspr,fname_ssmt,
     &  fname_sweobs,fname_sweobs_barnes_mask,iday_init,
     &  imonth_init,iyear_init,nyear,nobs_total,depth_assim,
     &  fname_sden)

      implicit none

      include 'snowmodel.inc'

      real deltax,deltay,undef,dn,beta,areas_flag,swe_count,
     &  cf_min,ro_avg
      real sprec_ratio(max_obs_dates),smelt_ratio(max_obs_dates)
      real corr_factor(nx_max,ny_max,max_obs_dates)
      real swe_tmp(nx_max,ny_max),sum_sprec_tmp1(nx_max,ny_max),
     &  sum_sprec_tmp2(nx_max,ny_max),grid(nx_max,ny_max),
     &  areas_mask(nx_max,ny_max),sum_smelt_tmp1(nx_max,ny_max),
     &  sum_smelt_tmp2(nx_max,ny_max),ro_tmp(nx_max,ny_max)
      real corr_factor_tmp(nx_max*ny_max),swe_obs(nx_max*ny_max),
     &  swe_model(nx_max*ny_max),sumsprec_model(nx_max*ny_max),
     &  delta_old(nx_max*ny_max),obsid(nx_max*ny_max),
     &  sumsmelt_model(nx_max*ny_max),obsid_old(nx_max*ny_max),
     &  delta_old_tmp(nx_max*ny_max),depth_obs(nx_max*ny_max),
     &  ro_model(nx_max*ny_max)
c     real corr_offset(nx_max*ny_max)

      double precision xmn,ymn
      double precision xstn(nx_max*ny_max),ystn(nx_max*ny_max)

      integer iobs_rec(max_obs_dates)
      integer ii(nx_max*ny_max),jj(nx_max*ny_max)
      integer iobs_num,irec1,irec2,nobs_dates,nx,ny,i,j,ifill,
     &  iobsint,k,nstns,nstns_old,kk,local_assim_flag,iiyr,iimo,
     &  iidy,iobs_rec_tmp,iday_init,imonth_init,iyear_init,nyear,
     &  krec,nobs_total
      integer depth_assim,count

      character*80 fname_swed,fname_sspr,fname_ssmt,fname_sden
      character*80 fname_sweobs
      character*80 fname_sweobs_barnes_mask

c Perform some initialization steps.
      if (nyear.eq.1) then

c Define some of the constants and parameters used in the data
c   assimilation.  ifill should be = 1; in that case undef is not
c   really used (so it does not have to be the same as defined
c   in the .par file).
        undef = -9999.0
        ifill = 1
        iobsint = 0

c Open the swe input file.
        open (unit=61,file=fname_sweobs)

c Open a file to write some basic correction factor information.
c   This just saves information that the user might want to look
c   at.
        open (unit=77,file='data/corr_factor.txt')

c Open an output file for the correction factor array.
        open(62,file='data/corr_factor.gdat',
     &    form='unformatted',access='direct',recl=4*nx*ny)

c Initialize the number-of-observations counter.
        nobs_total = 0

      endif

c Read the number of observation dates for this year.
      read(61,*) nobs_dates

c If you have observations for this year, generate the correction
c   factors.  For the case of no observations for this year, set
c   the correction factor equal to 1.0.
      if (nobs_dates.gt.0) then

c Loop through the observation dates.
        do iobs_num=1,nobs_dates

c Increment the number-of-observations counter.
          nobs_total = nobs_total + 1

c Read the date corresponding to this observation.
          read(61,*) iiyr,iimo,iidy

c Convert this date to the corresponding record number in the
c   original SnowModel output files.  Note that this has assumed
c   that daily data files were written out since the start of
c   the simulation.
          call get_obs_record(iday_init,imonth_init,iyear_init,
     &      iidy,iimo,iiyr,iobs_rec_tmp)
          iobs_rec(iobs_num) = iobs_rec_tmp

c For this observation date, read in the data describing the
c   location and swe values for each observation.  For areas
c   simulations, xstn, and ystn correspond to the center of the
c   area domain and they are not really used.
          read(61,*) nstns
          do k=1,nstns
            if (depth_assim.eq.1) then
              read(61,*) obsid(k),xstn(k),ystn(k),depth_obs(k)
            else
              read(61,*) obsid(k),xstn(k),ystn(k),swe_obs(k)
            endif
          enddo

c Convert the x and y locations to (ii,jj) locations.
          do k=1,nstns
            ii(k) = 1 + nint((xstn(k) - xmn) / deltax)
            jj(k) = 1 + nint((ystn(k) - ymn) / deltay)
          enddo


c If you do a data assimilation run from start to finish, it is
c   not required to close and reopen these files.  But if you are
c   doing a history restart then these files are no longer open
c   so you must do this.
          close (238)
          close (239)
          close (240)
          close (237)

c Open the required inputs from the initial assimilation loop.
c   Open swe depth (swe_depth).
c     /outputs/wo_assim/swed.gdat is unit 238 in outputs_user.f
c   Open sum snow precip (sum_sprec).
c     /outputs/wo_assim/sspr.gdat is unit 239 in outputs_user.f
c   Open sum snow melt (sum_smelt).
c     /outputs/wo_assim/ssmt.gdat is unit 240 in outputs_user.f
          open (238,file=fname_swed,
     &      form='unformatted',access='direct',recl=4*1*nx*ny)
          open (239,file=fname_sspr,
     &      form='unformatted',access='direct',recl=4*1*nx*ny)
          open (240,file=fname_ssmt,
     &      form='unformatted',access='direct',recl=4*1*nx*ny)
          print *,fname_ssmt
          print *,fname_sden
          open (237,file=fname_sden,
     &      form='unformatted',access='direct',recl=4*1*nx*ny)

c Read the model output for the first observation time.
          if (iobs_num.eq.1) then
            irec1 = iobs_rec(iobs_num)
            read(238,rec=irec1) ((swe_tmp(i,j),i=1,nx),j=1,ny)
            read(239,rec=irec1) ((sum_sprec_tmp1(i,j),i=1,nx),j=1,ny)
            read(240,rec=irec1) ((sum_smelt_tmp1(i,j),i=1,nx),j=1,ny)
            read(237,rec=irec1) ((ro_tmp(i,j),i=1,nx),j=1,ny)

            count = 0
            ro_avg = 0
            do i=1,nx
              do j=1,ny
                if (ro_tmp(i,j).ge.0) then
                  ro_avg = ro_avg + ro_tmp(i,j)
                  count = count + 1
                endif
              enddo
            enddo
            if (count.gt.0) ro_avg = ro_avg/count
              
c            print *,swe_tmp(1,1), sum_sprec_tmp1(1,1),
c     &        sum_smelt_tmp1(1,1)
c            stop

c For points, just pull the data at the appropriate grid cell.
c   For areas, average the data over the masked out area for each
c   'station'.
            do k=1,nstns
              if (areas_flag.eq.0.0) then
                swe_model(k) = swe_tmp(ii(k),jj(k))
                ro_model(k) = ro_tmp(ii(k),jj(k))
                sumsprec_model(k) = sum_sprec_tmp1(ii(k),jj(k))
                sumsmelt_model(k) = sum_smelt_tmp1(ii(k),jj(k))
              elseif (areas_flag.eq.1.0) then
                swe_model(k) = 0.0
                ro_model(k) = 0.0
                sumsprec_model(k) = 0.0
                sumsmelt_model(k) = 0.0
                swe_count = 0.0
                do j=1,ny
                  do i=1,nx
                    if (areas_mask(i,j).eq.obsid(k)) then
c The following is used if the mask changes with observation time.
c                   if (areas_mask(i,j,iobs_num).eq.obsid(k)) then
                      swe_count = swe_count + 1.0
                      swe_model(k) = swe_model(k) + swe_tmp(i,j)
                      ro_model(k) = ro_model(k) + ro_tmp(i,j)
                      sumsprec_model(k) = sumsprec_model(k) +
     &                  sum_sprec_tmp1(i,j)
                      sumsmelt_model(k) = sumsmelt_model(k) +
     &                  sum_smelt_tmp1(i,j)
                    endif
                  enddo
                enddo
                swe_model(k) = swe_model(k) / swe_count
                ro_model(k) = ro_model(k) / swe_count
                sumsprec_model(k) = sumsprec_model(k) / swe_count
                sumsmelt_model(k) = sumsmelt_model(k) / swe_count
              endif
            enddo
          endif

c Read the model output for any additional observation times (irec1
c   = current obs time, irec2 = previous obs time).
          if (iobs_num.gt.1) then
            irec1 = iobs_rec(iobs_num)
            irec2 = iobs_rec(iobs_num-1)
            read(238,rec=irec1) ((swe_tmp(i,j),i=1,nx),j=1,ny)
            read(239,rec=irec1) ((sum_sprec_tmp1(i,j),i=1,nx),j=1,ny)
            read(239,rec=irec2) ((sum_sprec_tmp2(i,j),i=1,nx),j=1,ny)
            read(240,rec=irec1) ((sum_smelt_tmp1(i,j),i=1,nx),j=1,ny)
            read(240,rec=irec2) ((sum_smelt_tmp2(i,j),i=1,nx),j=1,ny)
            read(237,rec=irec1) ((ro_tmp(i,j),i=1,nx),j=1,ny)

            count = 0
            ro_avg = 0
            do i=1,nx
              do j=1,ny
                if (ro_tmp(i,j).ge.0) then
                  ro_avg = ro_avg + ro_tmp(i,j)
                  count = count + 1
                endif
              enddo
            enddo
            if (count.gt.0) ro_avg = ro_avg/count

c For points, just pull the data at the appropriate grid cell.
c   For areas, average the data over the masked out area for each
c   'station'.
            do k=1,nstns
              if (areas_flag.eq.0.0) then
                swe_model(k) = swe_tmp(ii(k),jj(k))
                ro_model(k) = ro_tmp(ii(k),jj(k))
                sumsprec_model(k) = sum_sprec_tmp1(ii(k),jj(k)) -
     &            sum_sprec_tmp2(ii(k),jj(k))
                sumsmelt_model(k) = sum_smelt_tmp1(ii(k),jj(k)) -
     &            sum_smelt_tmp2(ii(k),jj(k))
              elseif (areas_flag.eq.1.0) then
                swe_model(k) = 0.0
                ro_model(k) = 0.0
                sumsprec_model(k) = 0.0
                sumsmelt_model(k) = 0.0
                swe_count = 0.0
                do j=1,ny
                  do i=1,nx
                    if (areas_mask(i,j).eq.obsid(k)) then
c The following is used if the mask changes with observation time.
c                   if (areas_mask(i,j,iobs_num).eq.obsid(k)) then
                      swe_count = swe_count + 1.0
                      swe_model(k) = swe_model(k) + swe_tmp(i,j)
                      ro_model(k) = ro_model(k) + ro_tmp(i,j)
                      sumsprec_model(k) = sumsprec_model(k) +
     &                  sum_sprec_tmp1(i,j) - sum_sprec_tmp2(i,j)
                      sumsmelt_model(k) = sumsmelt_model(k) +
     &                  sum_smelt_tmp1(i,j) - sum_smelt_tmp2(i,j)
                    endif
                  enddo
                enddo
                swe_model(k) = swe_model(k) / swe_count
                ro_model(k) = ro_model(k) / swe_count
                sumsprec_model(k) = sumsprec_model(k) / swe_count
                sumsmelt_model(k) = sumsmelt_model(k) / swe_count
              endif
            enddo
          endif

c To avoid a divide by zero later on, make sure sumsprec_model and
c   sumsmelt_model are not both zero.
          do k=1,nstns
            sumsprec_model(k) = sumsprec_model(k) + 1.0e-6
          enddo

c Determine whether we will adjust the precipitation or melt.  To
c   do this, calculate the relative contributions of precipitation
c   and melt inputs for this correction period.  This can be
c   different for each observation interval.  Calculate the average
c   over all of the stations/areas in the domain.
          sprec_ratio(iobs_num) = 0.0
          smelt_ratio(iobs_num) = 0.0
          do k=1,nstns
            sprec_ratio(iobs_num) = sprec_ratio(iobs_num) +
     &        sumsprec_model(k) / (sumsprec_model(k)+sumsmelt_model(k))
            smelt_ratio(iobs_num) = smelt_ratio(iobs_num) +
     &        sumsmelt_model(k) / (sumsprec_model(k)+sumsmelt_model(k))
          enddo
          sprec_ratio(iobs_num) = sprec_ratio(iobs_num) / real(nstns)
          smelt_ratio(iobs_num) = smelt_ratio(iobs_num) / real(nstns)

c Initialize the delta swe variable.
          if (iobs_num.eq.1) then
            do k=1,nstns
              delta_old(k) = 0.0
            enddo
          else
            do k=1,nstns
              delta_old(k) = 0.0
            enddo
            do k=1,nstns
              do kk=1,nstns_old
                if(obsid(k).eq.obsid_old(kk))
     &            delta_old(k) = delta_old_tmp(kk)
              enddo
            enddo
c           write (77,*)
c           do k=1,nstns
c             write (77,*) 'k, delta_old(k)',k,100.*delta_old(k)
c           enddo
c           write (77,*)
          endif

c Calculate the correction factor to be used in the next model
c   iteration.  Let the correction factor equal 1.0 during
c   periods where we have no swe observations.  Also, note that the
c   reason for the delta_old variable is to account for the fact
c   that that delta will be fixed with the previous date correction
c   time period.  This is one of the things that allows the
c   correction to be done in two model iterations.
c If sumsprec_model or sumsmelt_model are too small to be used in
c   the assimilation (like less than 1 mm), set corr_factor_tmp = 1.0
c   so no adjustments are performed for this observation interval.
          cf_min = 0.1
          do k=1,nstns
            if (depth_assim.eq.1) then
              if (ro_model(k).le.0) then
                swe_obs(k) = ro_avg*depth_obs(k)/1000.0
              else
                swe_obs(k) = ro_model(k)*depth_obs(k)/1000.0
              endif
            endif
            if (sprec_ratio(iobs_num).ge.smelt_ratio(iobs_num)) then
              if (sumsprec_model(k).lt.1.0e-3) then
                corr_factor_tmp(k) = 1.0
              else
                corr_factor_tmp(k) = 1.0 +
     &            (swe_obs(k) - swe_model(k) - delta_old(k)) /
     &            sumsprec_model(k)
                corr_factor_tmp(k) = max(cf_min,corr_factor_tmp(k))
              endif
            else
              if (sumsmelt_model(k).lt.1.0e-3) then
                corr_factor_tmp(k) = 1.0
              else
                corr_factor_tmp(k) = 1.0 +
     &            (swe_model(k) - swe_obs(k) + delta_old(k)) /
     &            sumsmelt_model(k)
                corr_factor_tmp(k) = max(cf_min,corr_factor_tmp(k))
              endif
            endif
c Save some information about the model calculations.
c           write (77,*) '---'
c           write (77,*) k,swe_obs(k)
c           write (77,*) k,swe_model(k)
c           write (77,*) k,delta_old(k)
c           write (77,*) k,swe_obs(k)-swe_model(k)-delta_old(k)
c           write (77,*) k,sumsprec_model(k)
c           write (77,*) k,sumsmelt_model(k)
c           write (77,*) k,corr_factor_tmp(k)
c           write (77,*) '---'

c Save some data from this observation time for use at the next
c   observation time.
            nstns_old = nstns
            obsid_old(k) = obsid(k)
            delta_old_tmp(k) = swe_obs(k) - swe_model(k)
          enddo

c Now that I have the correction factors calculated at each
c   observation point, interpolate those over the simulation domain.

c Use the barnes oi scheme to create the distribution. If there is
c   only a single station, distribute those data uniformly over
c   the domain.
          if (nstns.ge.2) then
            call get_dn(nx,ny,deltax,deltay,nstns,dn,iobsint)

c Modify the size of dn.
            dn = beta * dn

            call barnes_oi(nx,ny,deltax,deltay,xmn,ymn,
     &        nstns,xstn,ystn,corr_factor_tmp,dn,grid,undef,ifill)
          elseif (nstns.eq.1) then
            call single_stn(nx,ny,nstns,corr_factor_tmp,grid)
          endif

c The following calculations are done if you want to implement the
c   special case where you limit the data assimilation corrections
c   to a certain area of your simulation domain.  Edits to this
c   subroutine will certainly be required to make it work for your
c   specific application.  This subroutine generates correction
c   factors with 1.0's in places outside the local obs influences.
          if (local_assim_flag.eq.1) then
            call mk_local_cfs(nx,ny,undef,xmn,ymn,deltax,deltay,
     &        fname_sweobs,fname_sweobs_barnes_mask,nobs_dates,
     &        corr_factor_tmp,beta,iobsint,ifill)
          endif

c Define the correction surface record that corresponds to this
c   year and observation.
          krec = nobs_total + (nyear - 1)
          if (krec.gt.max_obs_dates) then
            print *, 'max_obs_dates must be increased in snowmodel.inc'
            print *, 'krec = ',krec,'  max_obs_dates = ',max_obs_dates
            stop
          endif

c Use the gridded output file to build the corr_factor array.
          do j=1,ny
            do i=1,nx
              corr_factor(i,j,krec) = grid(i,j)
              corr_factor(i,j,krec) =
     &          max(cf_min,corr_factor(i,j,krec))
            enddo
          enddo

c Note that the interpolation scheme may have produced correction
c   factors that do not produce exact matches with the
c   observations (like happens with the case of having a single
c   data point).  If you are interested, calculate the difference
c   between the exact value and the actual calculated value, and
c   then write it out as done below.
c         do k=1,nstns
c           if (sprec_ratio(iobs_num).ge.smelt_ratio(iobs_num)) then
c             corr_offset(k) = sumsprec_model(k) *
c    &          (corr_factor(ii(k),jj(k),iobs_num) - corr_factor_tmp(k))
c           else
c             corr_offset(k) = sumsmelt_model(k) *
c    &          (corr_factor(ii(k),jj(k),iobs_num) - corr_factor_tmp(k))
c           endif
c         enddo

c Write some information to the text file.
          write (77,*) '***************************************'
          write (77,*) ' sprec_ratio =',sprec_ratio(iobs_num),
     &      '  smelt_ratio =',smelt_ratio(iobs_num)
          write (77,*)
          do k=1,nstns
            write (77,*) k,' swe diff =',
     &        100.0*abs(swe_obs(k)-swe_model(k)),' SWE OBS =',
     &        100.0*swe_obs(k)
            write (77,*) 'sumsprec =',sumsprec_model(k)*100.,
     &        '  SWE MODEL =',swe_model(k)*100.
            write (77,*) 'iobs_num =',iobs_num,
     &        '  CORR_FACTOR =',corr_factor_tmp(k)
cc          write (77,*) 'corr_offset =',100.*corr_offset(k),
cc   &        '  ij',ii(k),jj(k)
cc          write (77,*) '     delta_old =',100.*delta_old(k),
cc   &        '      corr fact used =',corr_factor(ii(k),jj(k),iobs_num)
c           write (77,*)
c           write (77,*) k,' sumsprec_model(k) =',sumsprec_model(k)
c           write (77,*) k,' sumsmelt_model(k) =',sumsmelt_model(k)
c           write (77,*)
          enddo
          write (77,*) '***************************************'

c Write the output data to a grads file.
          write(62,rec=krec) ((corr_factor(i,j,krec),i=1,nx),j=1,ny)

        enddo

c Fill corr_factor with 1.0 for the period following the last obs
c   date in the current year.  This is also required for the history
c   restart to work correctly.  Without the history restart this was
c   already done as part of the model initialization.
        if (krec+1.gt.max_obs_dates) then
          print *, 'max_obs_dates must be increased in snowmodel.inc'
          print *, 'krec+1 = ',krec+1,'  max_obs_dates = ',max_obs_dates
          stop
        endif
        do j=1,ny
          do i=1,nx
            corr_factor(i,j,krec+1) = 1.0
          enddo
        enddo

        write(62,rec=krec+1) ((corr_factor(i,j,krec+1),i=1,nx),j=1,ny)

c The met, topo, and veg files must be closed for the next model
c   iteration.
        close (20)
        close (37)
        close (38)

        close (238)
        close (239)
        close (240)
        close (237)

      else

c For the case of no observations for this year, set the correction
c   factor equal to 1.0.
        krec = nobs_total + nyear
        if (krec.gt.max_obs_dates) then
          print *, 'max_obs_dates must be increased in snowmodel.inc'
          print *, 'krec = ',krec,'  max_obs_dates = ',max_obs_dates
          stop
        endif

        do j=1,ny
          do i=1,nx
            corr_factor(i,j,krec) = 1.0
          enddo
        enddo

        write(62,rec=krec) ((corr_factor(i,j,krec),i=1,nx),j=1,ny)

      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine corr_factor_index(nobs_dates,icorr_factor_index,
     &  iobs_rec,max_iter,sprec_ratio,smelt_ratio,print_inc,
     &  iday_init,imonth_init,iyear_init,nyear,nobs_total_cfi,
     &  nyears,ihrestart_flag)

      implicit none

      include 'snowmodel.inc'

      integer icorr_factor_index(max_time_steps)

      integer kk,istart,iend,nobs_dates,iter,max_iter,krec,nyear,
     &  nobs_total_cfi,ioptn,julian_start,iiyr,julian_end,iday_init,
     &  imonth_init,iyear_init,nyears,ihrestart_flag
      integer iobs_rec(max_obs_dates)
      real sprec_ratio(max_obs_dates),smelt_ratio(max_obs_dates) 
      real print_inc

      if (ihrestart_flag.eq.0) print_inc = 24.0

c Initialize the number-of-observations counter.
      if (nyear.eq.1) nobs_total_cfi = 0

c Build an array indicating the appropriate correction factor to
c   use at each time step during the simulation.
      if (nobs_dates.gt.0) then

c Loop through the observation dates.
        do kk=1,nobs_dates+1

c Increment the number-of-observations counter.
          if (kk.le.nobs_dates) nobs_total_cfi = nobs_total_cfi + 1

c FIRST, FROM THE SIMULATION START UNTIL THE FIRST OBSERVATION, FOR
c   EACH YEAR.
          if (kk.eq.1) then

c Here istart equals the first model time step of each year.
            ioptn = 3
            call calndr (ioptn,iday_init,imonth_init,iyear_init,
     &        julian_start)

c Find the Julian day for this data record.
            iiyr = iyear_init + (nyear - 1)
            call calndr (ioptn,iday_init,imonth_init,iiyr,julian_end)

c Calculate istart.
            istart = 1 + (julian_end - julian_start) * nint(print_inc)

c Here iend equals the model time step corresponding to the end of
c   the first observation day.  Take advantage of the data-file
c   record that was calculated before.  Convert from daily data
c   output records to model time steps using print_inc.
            iend = iobs_rec(kk) * nint(print_inc)

c Define the corr_factor data array record.
            krec = nobs_total_cfi + (nyear - 1)
c Fill the index for each model time step.
            do iter=istart,iend
              if (sprec_ratio(kk).ge.smelt_ratio(kk)) then
                icorr_factor_index(iter) = krec
              else
                icorr_factor_index(iter) = -krec
              endif
            enddo

c SECOND, BETWEEN THE LAST OBSERVATION AND THE END OF THE SIMULATION.
          elseif (kk.eq.nobs_dates+1) then
            istart = iobs_rec(kk-1) * nint(print_inc) + 1

c Here iend equals the last time step in the year of interest.
c Find the Julian day for this data record.
            iiyr = iyear_init + nyear
            call calndr (ioptn,iday_init,imonth_init,iiyr,julian_end)

c Calculate iend for this year.
            iend = (julian_end - julian_start) * nint(print_inc)
            iend = min(iend,max_iter)

c Define the corr_factor data array record.
            krec = nobs_total_cfi + nyear

c Fill the index for each model time step.
            do iter=istart,iend
              icorr_factor_index(iter) = krec
            enddo

c THIRD, ANY PERIODS BETWEEN OBSERVATIONS.
          else
            istart = iobs_rec(kk-1) * nint(print_inc) + 1
            iend = iobs_rec(kk) * nint(print_inc)

c Define the corr_factor data array record.
            krec = nobs_total_cfi + (nyear - 1)

c Fill the index for each model time step.
            do iter=istart,iend
              if (sprec_ratio(kk).ge.smelt_ratio(kk)) then
                icorr_factor_index(iter) = krec
              else
                icorr_factor_index(iter) = -krec
              endif
            enddo
          endif
        enddo

      else

c Create an array indes for the case of no observations for this
c   year.  Here istart equals the first model time step of each year.
        ioptn = 3
        call calndr (ioptn,iday_init,imonth_init,iyear_init,
     &    julian_start)

c Find the Julian day for this data record.
        iiyr = iyear_init + (nyear - 1)
        call calndr (ioptn,iday_init,imonth_init,iiyr,julian_end)

c Calculate istart.
        istart = 1 + (julian_end - julian_start) * nint(print_inc)

c Calculate iend for this year.
        iiyr = iyear_init + nyear
        call calndr (ioptn,iday_init,imonth_init,iiyr,julian_end)
        iend = (julian_end - julian_start) * nint(print_inc)
        iend = min(iend,max_iter)

c Define the corr_factor data array record.
        krec = nobs_total_cfi + nyear

c Fill the index for each model time step.
        do iter=istart,iend
          icorr_factor_index(iter) = krec
        enddo

      endif

c SAVE A VERSION OF THE INDEX THAT THE USER CAN EASILY LOOK AT.
      if (nyear.eq.nyears) then
        print *
        do iter=1,max_iter
          write (77,*) iter,icorr_factor_index(iter)
        enddo
        print *
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mk_local_cfs(nx,ny,undef,xmn,ymn,deltax,deltay,
     &  fname_sweobs,fname_sweobs_barnes_mask,nobs_dates,
     &  corr_factor_tmp,beta,iobsint,ifill)

      implicit none

      include 'snowmodel.inc'

      integer nstns,k,nx,ny,i,j,ntstations,icount,nobs_dates,ifill,
     &  iobsint,nstns2

      real dummy,stnid,undef,deltax,deltay,beta,dn
      real xmask(nx_max,ny_max),grid(nx_max,ny_max)
      real corr_factor_tmp(nx_max*ny_max),cf_tmp2(nx_max*ny_max),
     &  obsid(nx_max*ny_max)

      real var(nstns_max)

      double precision xmn,ymn
      double precision yg(nx_max,ny_max),xg(nx_max,ny_max)
      double precision xstn(nstns_max),ystn(nstns_max)

      character*80 fname_sweobs
      character*80 fname_sweobs_barnes_mask

      print *
      print *,'You are doing local assimilation, this requires'
      print *,'an observation mask to have been generated before'
      print *,'the model run and saved in the file called:'
      print *,'fname_sweobs_barnes_mask'
      print *
      print *,'This was also not tested after some big changes to'
      print *,'the data assimilation code.  So I strongly suggest'
      print *,'you make sure this is doing what you want when you'
      print *,'first start using it.'
      print *
      stop
      if (nobs_dates.gt.1) then
        print *,'THIS HAS NOT BEEN MADE TO WORK WITH MORE THAN'
        print *,'ONE OBS TIME.'
        stop
      endif
      print *

c Save the correction factors for the local-influence assimilation
c   scheme.
      open (unit=78,file='data/cf.txt')
      do k=1,nstns
        write (78,*) corr_factor_tmp(k)
      enddo
      close (78)

c These are the obs and correction factors calculated from the
c   first loop in SnowModel.
      open (721,file=fname_sweobs)
      open (731,file='data/cf.txt')

      read (721,*) nstns
      do k=1,nstns
        read (721,*) stnid,xstn(k),ystn(k),dummy
        read (731,*) var(k)
c       if (var(k).gt.1.5) then
c         print *, 'cf > 1.5 found; setting to 1.5',k,var(k)
c         var(k) = 1.5
c       endif
c       if (var(k).lt.0.5) then
c         print *, 'cf < 0.5 found; setting to 0.5',k,var(k)
c         var(k) = 0.5
c       endif
c       var(k) = min(1.5,var(k))
c       var(k) = max(0.5,var(k))
      enddo
      close (721)
      close (731)

c Create a collection of 'stations' with correction factors of
c   1.0 in areas outside of our traverse regions.
      open(741,file=fname_sweobs_barnes_mask,
     &  form='unformatted',access='direct',recl=4*nx*ny)
      read (741,rec=1) ((xmask(i,j),i=1,nx),j=1,ny)
      close (741)

c Create an array of e, n coordinates for this domain.
      do j=1,ny
        do i=1,nx
          xg(i,j) = xmn + deltax * (real(i) - 1.0)
          yg(i,j) = ymn + deltay * (real(j) - 1.0)
        enddo
      enddo

      do j=1,ny
        do i=1,nx
          if (xmask(i,j).ne.undef) then
            xg(i,j) = undef
            yg(i,j) = undef
          endif
        enddo
      enddo

c Count how many cf=1.0 'stations' you are going to end up with.
      icount = 0
      do j=1,ny,100
        do i=1,nx,100
          if (xg(i,j).ne.undef) then
            icount = icount + 1
          endif
        enddo
      enddo

c Write out the original stations.
      open (761,file='data/cf_with_mask.txt')
      ntstations = nstns + icount
      write (761,88) ntstations
      do k=1,nstns
        write (761,89) k,xstn(k),ystn(k),var(k)
      enddo

c Write out the cf=1.0 stations.
      icount = 0
      do j=1,ny,100
        do i=1,nx,100
          if (xg(i,j).ne.undef) then
            icount = icount + 1
            write (761,89) icount+1000,xg(i,j),yg(i,j),1.0
          endif
        enddo
      enddo
      close (761)

c Read in the new local cf data.
      open (79,file='data/cf_with_mask.txt')
      read (79,*) nstns2
      do k=1,nstns2
        read (79,*) obsid(k),xstn(k),ystn(k),cf_tmp2(k)
      enddo
      close (79)

      call get_dn(nx,ny,deltax,deltay,nstns2,dn,iobsint)

      dn = beta * dn

      call barnes_oi(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns2,xstn,ystn,cf_tmp2,dn,grid,undef,ifill)

c Write the output data to a grads file.
      open(511,file='data/corr_factor_w-mask.gdat',
     &  form='unformatted',access='direct',recl=4*nx*ny)
      write(511,rec=1) ((grid(i,j),i=1,nx),j=1,ny)
      close (511)

  88  format (i10)
  89  format (i10,2f15.1,f10.4)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_obs_record(iday_init,imonth_init,iyear_init,
     &  iidy,iimo,iiyr,iobs_rec_tmp)

      implicit none

      integer ioptn,iday_init,imonth_init,iyear_init,julian_start,
     &  iidy,iimo,iiyr,julian_end,iobs_rec_tmp

c Find the Julian day at the start of the model run.
      ioptn = 3
      call calndr (ioptn,iday_init,imonth_init,iyear_init,julian_start)

c Find the Julian day for this data record.
      call calndr (ioptn,iidy,iimo,iiyr,julian_end)

c Calculate the day of simulation for this data record.  This is the
c   same as the output file data record.
      iobs_rec_tmp = julian_end - julian_start + 1

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c file = calendar.f  version 1.0
c
c Note from Glen: This is the version that should be used as
c   part of my computer programs!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c This program performs various date conversion calculations
c using subroutine calndr().
c On a 32-bit computer, it can handle any date between roughly
c 5 million BC and 5 million AD.  This limitation is due to
c the range of integers that can be expressed with 32 bits.
c The algorithm has no limitation.
c
c Using function idaywk(), the day of the week is computed
c along with the answer to the user's calendar calculation.
c
c External routines called:
c calndr  calendar conversions
c idaywk  day of the week determination
c
c Portability
c This routine is coded to Fortran 77 standards except that
c lower case is used.
c
c Copyright (C) 1999 Jon Ahlquist.
c Issued under the second GNU General Public License.
c See www.gnu.org for details.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c If you find any errors, please notify:
c Jon Ahlquist <ahlquist@met.fsu.edu>   
c Dept of Meteorology
c Florida State University
c Tallahassee, FL 32306-4520
c 15 March 1999.
c
c----------
c
c Declare variables.
c     implicit     none
c     integer      day, daynum, ioptn, julian, month, year
c     character*9  daynam(0:6)
c
c Declare the integer function used to compute the day of the week.
c     integer  idaywk
c
c Define the day names.
c     data  daynam /'Sunday',   'Monday', 'Tuesday', 'Wednesday',
c    &              'Thursday', 'Friday', 'Saturday'/
c
c Variables and their meanings
c day     day of the month.
c daynam  array of day names.  (daynam(0)='Sunday', daynam(1)='Monday',
c            ..., daynam(6)='Saturday')
c daynum  day number during the year.  (1 for 1 January, 2 for
c            2 January, 32 for 1 February, etc.)
c idaywk  integer function that returns an integer counter indicating
c            the day of the week, where 0 refers to Sunday, 1 to Monday,
c            up to 6 for Saturday.
c ioptn   option indicator where 0 < abs(ioptn) < 6.
c            See below and especially subroutine calndr for details.
c julian  Julian Day number.
c month   month counter (1=January, 2=February, ..., 12=December)
c year    year expressed with ALL digits.  DO NOT abbreviate years
c            by using only the last two digits.
c
c----------
c
      subroutine calndr (ioptn, iday, month, iyear, idayct)
c
c----------
c
c CALNDR = CALeNDaR conversions, version 1.0
c
c Input variable specifying the desired calendar conversion option.
      integer ioptn
c
c Input/Output variables (sometimes input, sometimes output,
c depending on the value of the desired option, ioptn.)
      integer  iday, month, iyear, idayct
c
c----------
c
c Subroutine calndr() performs calendar calculations using either
c the standard Gregorian calendar or the old Julian calendar.
c This subroutine extends the definitions of these calendar systems
c to any arbitrary year.  The algorithms in this subroutine
c will work with any date in the past or future,
c but overflows will occur if the numbers are sufficiently large.
c For a computer using a 32-bit integer, this routine can handle
c any date between roughly 5.8 million BC and 5.8 million AD
c without experiencing overflow during calculations.
c
c No external functions or subroutines are called.
c
c----------
c
c INPUT/OUTPUT ARGUMENTS FOR SUBROUTINE CALNDR()
c
c "ioptn" is the desired calendar conversion option explained below.
c Positive option values use the standard modern Gregorian calendar.
c Negative option values use the old Julian calendar which was the
c standard in Europe from its institution by Julius Caesar in 45 BC
c until at least 4 October 1582.  The Gregorian and Julian calendars
c are explained further below.
c
c (iday,month,iyear) is a calendar date where "iday" is the day of
c the month, "month" is 1 for January, 2 for February, etc.,
c and "iyear" is the year.  If the year is 1968 AD, enter iyear=1968,
c since iyear=68 would refer to 68 AD.
c For BC years, iyear should be negative, so 45 BC would be iyear=-45.
c By convention, there is no year 0 under the BC/AD year numbering
c scheme.  That is, years proceed as 2 BC, 1 BC, 1 AD, 2 AD, etc.,
c without including 0.  Subroutine calndr() will print an error message
c and stop if you specify iyear=0.
c
c "idayct" is a day count.  It is either the day number during the
c specified year or the Julian Day number, depending on the value
c of ioptn.  By day number during the specified year, we mean
c idayct=1 on 1 January, idayct=32 on 1 February, etc., to idayct=365
c or 366 on 31 December, depending on whether the specified year
c is a leap year.
c
c The values of input variables are not changed by this subroutine.
c
c
c ALLOWABLE VALUES FOR "IOPTN" and the conversions they invoke.
c Positive option values ( 1 to  5) use the standard Gregorian calendar.
c Negative option values (-1 to -5) use the old      Julian    calendar.
c
c Absolute
c  value
c of ioptn   Input variable(s)     Output variable(s)
c
c    1       iday,month,iyear      idayct
c Given a calendar date (iday,month,iyear), compute the day number
c (idayct) during the year, where 1 January is day number 1 and
c 31 December is day number 365 or 366, depending on whether it is
c a leap year.
c
c    2       idayct,iyear          iday,month
c Given the day number of the year (idayct) and the year (iyear),
c compute the day of the month (iday) and the month (month).
c
c    3       iday,month,iyear      idayct
c Given a calendar date (iday,month,iyear), compute the Julian Day
c number (idayct) that starts at noon of the calendar date specified.
c
c    4       idayct                iday,month,iyear
c Given the Julian Day number (idayct) that starts at noon,
c compute the corresponding calendar date (iday,month,iyear).
c
c    5       idayct                iday,month,iyear
c Given the Julian Day number (idayct) that starts at noon,
c compute the corresponding day number for the year (iday)
c and year (iyear).  On return from calndr(), "month" will always
c be set equal to 1 when ioptn=5.
c
c No inverse function is needed for ioptn=5 because it is
c available through option 3.  One simply calls calndr() with:
c ioptn = 3,
c iday  = day number of the year instead of day of the month,
c month = 1, and
c iyear = whatever the desired year is.
c
c----------
c
c EXAMPLES
c The first 6 examples are for the standard Gregorian calendar.
c All the examples deal with 15 October 1582, which was the first day
c of the Gregorian calendar.  15 October is the 288-th day of the year.
c Julian Day number 2299161 began at noon on 15 October 1582.
c
c Find the day number during the year on 15 October 1582
c     ioptn = 1
c     call calndr (ioptn, 15, 10, 1582,  idayct)
c calndr() should return idayct=288
c
c Find the day of the month and month for day 288 in year 1582.
c     ioptn = 2
c     call calndr (ioptn, iday, month, 1582, 288)
c calndr() should return iday=15 and month=10.
c
c Find the Julian Day number for 15 October 1582.
c     ioptn = 3
c     call calndr (ioptn, 15, 10, 1582, julian)
c calndr() should return julian=2299161
c
c Find the Julian Day number for day 288 during 1582 AD.
c When the input is day number of the year, one should specify month=1
c     ioptn = 3
c     call calndr (ioptn, 288, 1, 1582, julian)
c calndr() should return dayct=2299161
c
c Find the date for Julian Day number 2299161.
c     ioptn = 4
c     call calndr (ioptn, iday, month, iyear, 2299161)
c calndr() should return iday=15, month=10, and iyear=1582
c 
c Find the day number during the year (iday) and year
c for Julian Day number 2299161.
c     ioptn = 5
c     call calndr (ioptn, iday, month, iyear, 2299161)
c calndr() should return iday=288, month=1, iyear=1582
c
c Given 15 October 1582 under the Gregorian calendar,
c find the date (idayJ,imonthJ,iyearJ) under the Julian calendar.
c To do this, we call calndr() twice, using the Julian Day number
c as the intermediate value.
c     call calndr ( 3, 15,        10, 1582,    julian)
c     call calndr (-4, idayJ, monthJ, iyearJ,  julian)
c The first call to calndr() should return julian=2299161, and
c the second should return idayJ=5, monthJ=10, iyearJ=1582
c
c----------
c
c BASIC CALENDAR INFORMATION
c
c The Julian calendar was instituted by Julius Caesar in 45 BC.
c Every fourth year is a leap year in which February has 29 days.
c That is, the Julian calendar assumes that the year is exactly
c 365.25 days long.  Actually, the year is not quite this long.
c The modern Gregorian calendar remedies this by omitting leap years
c in years divisible by 100 except when the year is divisible by 400.
c Thus, 1700, 1800, and 1900 are leap years under the Julian calendar
c but not under the Gregorian calendar.  The years 1600 and 2000 are
c leap years under both the Julian and the Gregorian calendars.
c Other years divisible by 4 are leap years under both calendars,
c such as 1992, 1996, 2004, 2008, 2012, etc.  For BC years, we recall
c that year 0 was omitted, so 1 BC, 5 BC, 9 BC, 13 BC, etc., and 401 BC,
c 801 BC, 1201 BC, etc., are leap years under both calendars, while
c 101 BC, 201 BC, 301 BC, 501 BC, 601 BC, 701 BC, 901 BC, 1001 BC,
c 1101 BC, etc., are leap years under the Julian calendar but not
c the Gregorian calendar.
c
c The Gregorian calendar is named after Pope Gregory XIII.  He declared
c that the last day of the old Julian calendar would be Thursday,
c 4 October 1582 and that the following day, Friday, would be reckoned
c under the new calendar as 15 October 1582.  The jump of 10 days was
c included to make 21 March closer to the spring equinox.
c
c Only a few Catholic countries (Italy, Poland, Portugal, and Spain)
c switched to the Gregorian calendar on the day after 4 October 1582.
c It took other countries months to centuries to change to the
c Gregorian calendar.  For example, England's first day under the
c Gregorian calendar was 14 September 1752.  The same date applied to
c the entire British empire, including America.  Japan, Russia, and many
c eastern European countries did not change to the Gregorian calendar
c until the 20th century.  The last country to change was Turkey,
c which began using the Gregorian calendar on 1 January 1927.
c
c Therefore, between the years 1582 and 1926 AD, you must know
c the country in which an event was dated to interpret the date
c correctly.  In Sweden, there was even a year (1712) when February
c had 30 days.  Consult a book on calendars for more details
c about when various countries changed their calendars.
c
c DAY NUMBER DURING THE YEAR
c The day number during the year is simply a counter equal to 1 on
c 1 January, 32 on 1 February, etc., thorugh 365 or 366 on 31 December,
c depending on whether the year is a leap year.  Sometimes this is
c called the Julian Day, but that term is better reserved for the
c day counter explained below.
c
c JULIAN DAY NUMBER
c The Julian Day numbering system was designed by Joseph Scaliger
c in 1582 to remove ambiguity caused by varying calendar systems.
c The name "Julian Day" was chosen to honor Scaliger's father,
c Julius Caesar Scaliger (1484-1558), an Italian scholar and physician
c who lived in France.  Because Julian Day numbering was especially
c designed for astronomers, Julian Days begin at noon so that the day
c counter does not change in the middle of an astronmer's observing
c period.  Julian Day 0 began at noon on 1 January 4713 BC under the
c Julian calendar.  A modern reference point is that 23 May 1968
c (Gregorian calendar) was Julian Day 2,440,000.
c
c JULIAN DAY NUMBER EXAMPLES
c
c The table below shows a few Julian Day numbers and their corresponding
c dates, depending on which calendar is used.  A negative 'iyear' refers
c to BC (Before Christ).
c
c                     Julian Day under calendar:
c iday  month   iyear     Gregorian   Julian
c  24     11   -4714            0        -38
c   1      1   -4713           38          0
c   1      1       1      1721426    1721424
c   4     10    1582      2299150    2299160
c  15     10    1582      2299161    2299171
c   1      3    1600      2305508    2305518
c  23      5    1968      2440000    2440013
c   5      7    1998      2451000    2451013
c   1      3    2000      2451605    2451618
c   1      1    2001      2451911    2451924
c
c From this table, we can see that the 10 day difference between the
c two calendars in 1582 grew to 13 days by 1 March 1900, since 1900 was
c a leap year under the Julian calendar but not under the Gregorian
c calendar.  The gap will widen to 14 days after 1 March 2100 for the
c same reason.
c 
c----------
c
c PORTABILITY
c
c This subroutine is written in standard FORTRAN 77.
c It calls no external functions or subroutines and should run
c without problem on any computer having a 32-bit word or longer.
c 
c----------
c
c ALGORITHM
c
c The goal in coding calndr() was clear, clean code, not efficiency.
c Calendar calculations usually take a trivial fraction of the time
c in any program in which dates conversions are involved.
c Data analysis usually takes the most time.
c
c Standard algorithms are followed in this subroutine.  Internal to
c this subroutine, we use a year counter "jyear" such that
c  jyear=iyear   when iyear is positive
c       =iyear+1 when iyear is negative.
c Thus, jyear does not experience a 1 year jump like iyear does
c when going from BC to AD.  Specifically, jyear=0 when iyear=-1,
c i.e., when the year is 1 BC.
c
c For simplicity in dealing with February, inside this subroutine,
c we let the year begin on 1 March so that the adjustable month,
c February is the last month of the year.
c It is clear that the calendar used to work this way because the
c months September, October, November, and December refer to
c 7, 8, 9, and 10.  For consistency, jyear is incremented on 1 March
c rather than on 1 January.  Of course, everything is adjusted back to
c standard practice of years beginning on 1 January before answers
c are returned to the routine that calls calndr().
c
c Lastly, we use a trick to calculate the number of days from 1 March
c until the end of the month that precedes the specified month.
c That number of days is int(30.6001*(month+1))-122,
c where 30.6001 is used to avoid the possibility of round-off and
c truncation error.  For example, if 30.6 were used instead,
c 30.6*5 should be 153, but round-off error could make it 152.99999,
c which would then truncated to 152, causing an error of 1 day.
c
c Algorithm reference:
c Dershowitz, Nachum and Edward M. Reingold, 1990: Calendrical
c Calculations.  Software-Practice and Experience, vol. 20, number 9
c (September 1990), pp. 899-928.
c
c Copyright (C) 1999 Jon Ahlquist.
c Issued under the second GNU General Public License.
c See www.gnu.org for details.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c If you find any errors, please notify:
c Jon Ahlquist <ahlquist@met.fsu.edu>   
c Dept of Meteorology
c Florida State University
c Tallahassee, FL 32306-4520
c 15 March 1999.
c
c-----
c Declare internal variables.
      integer  jdref,  jmonth, jyear, leap,
     &         n1yr, n4yr, n100yr, n400yr,
     &         ndays, ndy400, ndy100, nyrs,
     &         yr400, yrref
c
c Explanation of all internal variables.
c jdref   Julian Day on which 1 March begins in the reference year.
c jmonth  Month counter which equals month+1 if month .gt. 2
c          or month+13 if month .le. 2.
c jyear   Year index,  jyear=iyear if iyear .gt. 0, jyear=iyear+1
c            if iyear .lt. 0.  Thus, jyear does not skip year 0
c            like iyear does between BC and AD years.
c leap    =1 if the year is a leap year, =0 if not.
c n1yr    Number of complete individual years between iyear and
c            the reference year after all 4, 100,
c            and 400 year periods have been removed.
c n4yr    Number of complete 4 year cycles between iyear and
c            the reference year after all 100 and 400 year periods
c            have been removed.
c n100yr  Number of complete 100 year periods between iyear and
c            the reference year after all 400 year periods
c            have been removed.
c n400yr  Number of complete 400 year periods between iyear and
c            the reference year.
c ndays   Number of days since 1 March during iyear.  (In intermediate
c            steps, it holds other day counts as well.)
c ndy400  Number of days in 400 years.  Under the Gregorian calendar,
c            this is 400*365 + 100 - 3 = 146097.  Under the Julian
c            calendar, this is 400*365 + 100 = 146100.
c ndy100  Number of days in 100 years,  Under the Gregorian calendar,
c            this is 100*365 + 24 = 36524.   Under the Julian calendar,
c            this is 100*365 + 25 = 36525.
c nyrs    Number of years from the beginning of yr400
c              to the beginning of jyear.  (Used for option +/-3).
c yr400   The largest multiple of 400 years that is .le. jyear.
c
c
c----------------------------------------------------------------
c Do preparation work.
c
c Look for out-of-range option values.
      if ((ioptn .eq. 0) .or. (abs(ioptn) .ge. 6)) then
         write(*,*)'For calndr(), you specified ioptn = ', ioptn
         write(*,*)
     &   'Allowable values are 1 to 5 for the Gregorian calendar'
         write(*,*)
     &   'and -1 to -5 for the Julian calendar.'
         stop
      endif
c
c Options 1-3 have "iyear" as an input value.
c Internally, we use variable "jyear" that does not have a jump
c from -1 (for 1 BC) to +1 (for 1 AD).
      if (abs(ioptn) .le. 3) then
         if (iyear .gt. 0) then
            jyear = iyear
         elseif (iyear .eq. 0) then
            write(*,*)
     &      'For calndr(), you specified the nonexistent year 0'
            stop
         else
            jyear = iyear + 1
         endif
c
c        Set "leap" equal to 0 if "jyear" is not a leap year
c        and equal to 1 if it is a leap year.
         leap = 0
         if ((jyear/4)*4 .eq. jyear) then
            leap = 1
         endif
         if ((ioptn .gt. 0)               .and.
     &       ((jyear/100)*100 .eq. jyear) .and.
     &       ((jyear/400)*400 .ne. jyear)      ) then
               leap = 0
         endif
      endif
c
c Options 3-5 involve Julian Day numbers, which need a reference year
c and the Julian Days that began at noon on 1 March of the reference
c year under the Gregorian and Julian calendars.  Any year for which
c "jyear" is divisible by 400 can be used as a reference year.
c We chose 1600 AD as the reference year because it is the closest
c multiple of 400 to the institution of the Gregorian calendar, making
c it relatively easy to compute the Julian Day for 1 March 1600
c given that, on 15 October 1582 under the Gregorian calendar,
c the Julian Day was 2299161.  Similarly, we need to do the same
c calculation for the Julian calendar.  We can compute this Julian
c Day knwoing that on 4 October 1582 under the Julian calendar,
c the Julian Day number was 2299160.  The details of these calculations
c is next. 
c    From 15 October until 1 March, the number of days is the remainder
c of October plus the days in November, December, January, and February:
c 17+30+31+31+28 = 137, so 1 March 1583 under the Gregorian calendar
c was Julian Day 2,299,298.  Because of the 10 day jump ahead at the
c switch from the Julian calendar to the Gregorian calendar, 1 March
c 1583 under the Julian calendar was Julian Day 2,299,308.  Making use
c of the rules for the two calendar systems, 1 March 1600 was Julian
c Day 2,299,298 + (1600-1583)*365 + 5 (due to leap years) =
c 2,305,508 under the Gregorian calendar and day 2,305,518 under the
c Julian calendar.
c    We also set the number of days in 400 years and 100 years.
c For reference, 400 years is 146097 days under the Gregorian calendar
c and 146100 days under the Julian calendar.  100 years is 36524 days
c under the Gregorian calendar and 36525 days under the Julian calendar.
      if (abs(ioptn) .ge. 3) then
c
c        Julian calendar values.
         yrref  =    1600
         jdref  = 2305518
c               = Julian Day reference value for the day that begins
c                 at noon on 1 March of the reference year "yrref".
         ndy400 = 400*365 + 100
         ndy100 = 100*365 +  25
c
c        Adjust for Gregorian calendar values.
         if (ioptn .gt. 0) then
            jdref  = jdref  - 10
            ndy400 = ndy400 -  3
            ndy100 = ndy100 -  1
         endif
      endif
c
c----------------------------------------------------------------
c OPTIONS -1 and +1:
c Given a calendar date (iday,month,iyear), compute the day number
c of the year (idayct), where 1 January is day number 1 and 31 December
c is day number 365 or 366, depending on whether it is a leap year.
      if (abs(ioptn) .eq. 1) then
c
c     Compute the day number during the year.
      if (month .le. 2) then
         idayct = iday + (month-1)*31
      else
         idayct = iday + int(30.6001 * (month+1)) - 63 + leap
      endif
c
c----------------------------------------------------------------
c OPTIONS -2 and +2:
c Given the day number of the year (idayct) and the year (iyear),
c compute the day of the month (iday) and the month (month).
      elseif (abs(ioptn) .eq. 2) then
c
      if (idayct .lt. 60+leap) then
         month  = (idayct-1)/31
         iday   = idayct - month*31
         month  = month + 1
      else
         ndays  = idayct - (60+leap)
c               = number of days past 1 March of the current year.
         jmonth = (10*(ndays+31))/306 + 3
c               = month counter, =4 for March, =5 for April, etc.
         iday   = (ndays+123) - int(30.6001*jmonth) 
         month  = jmonth - 1
      endif
c
c----------------------------------------------------------------
c OPTIONS -3 and +3:
c Given a calendar date (iday,month,iyear), compute the Julian Day
c number (idayct) that starts at noon.
      elseif (abs(ioptn) .eq. 3) then
c
c     Shift to a system where the year starts on 1 March, so January
c     and February belong to the preceding year.
c     Define jmonth=4 for March, =5 for April, ..., =15 for February.
      if (month .le. 2) then
        jyear  = jyear -  1
        jmonth = month + 13
      else
        jmonth = month +  1
      endif
c
c     Find the closest multiple of 400 years that is .le. jyear.
      yr400 = (jyear/400)*400
c           = multiple of 400 years at or less than jyear.
      if (jyear .lt. yr400) then
         yr400 = yr400 - 400
      endif
c
      n400yr = (yr400 - yrref)/400
c            = number of 400-year periods from yrref to yr400.
      nyrs   = jyear - yr400
c            = number of years from the beginning of yr400
c              to the beginning of jyear.
c
c     Compute the Julian Day number.
      idayct = iday + int(30.6001*jmonth) - 123 + 365*nyrs + nyrs/4
     &       + jdref + n400yr*ndy400
c
c     If we are using the Gregorian calendar, we must not count
c     every 100-th year as a leap year.  nyrs is less than 400 years,
c     so we do not need to consider the leap year that would occur if
c     nyrs were divisible by 400, i.e., we do not add nyrs/400.
      if (ioptn .gt. 0) then
         idayct = idayct - nyrs/100
      endif
c
c----------------------------------------------------------------
c OPTIONS -5, -4, +4, and +5:
c Given the Julian Day number (idayct) that starts at noon,
c compute the corresponding calendar date (iday,month,iyear)
c (abs(ioptn)=4) or day number during the year (abs(ioptn)=5).
      else
c
c     Create a new reference date which begins on the nearest
c     400-year cycle less than or equal to the Julian Day for 1 March
c     in the year in which the given Julian Day number (idayct) occurs.
      ndays  = idayct - jdref
      n400yr = ndays / ndy400
c            = integral number of 400-year periods separating
c              idayct and the reference date, jdref.
      jdref  = jdref + n400yr*ndy400
      if (jdref .gt. idayct) then
         n400yr = n400yr - 1
         jdref  = jdref  - ndy400
      endif
c
      ndays  = idayct - jdref
c            = number from the reference date to idayct.
c
      n100yr = min(ndays/ndy100, 3)
c            = number of complete 100-year periods
c              from the reference year to the current year.
c              The min() function is necessary to avoid n100yr=4
c              on 29 February of the last year in the 400-year cycle.
c
      ndays  = ndays - n100yr*ndy100
c            = remainder after removing an integral number of
c              100-year periods.
c
      n4yr   = ndays / 1461
c            = number of complete 4-year periods in the current century.
c              4 years consists of 4*365 + 1 = 1461 days.
c
      ndays  = ndays - n4yr*1461
c            = remainder after removing an integral number
c              of 4-year periods.
c
      n1yr   = min(ndays/365, 3)
c            = number of complete years since the last leap year.
c              The min() function is necessary to avoid n1yr=4
c              when the date is 29 February on a leap year,
c              in which case ndays=1460, and 1460/365 = 4.
c
      ndays  = ndays - 365*n1yr
c            = number of days so far in the current year,
c              where ndays=0 on 1 March.
c
      iyear  = n1yr + 4*n4yr + 100*n100yr + 400*n400yr + yrref 
c            = year, as counted in the standard way,
c              but relative to 1 March.
c
c At this point, we need to separate ioptn=abs(4), which seeks a
c calendar date, and ioptn=abs(5), which seeks the day number during
c the year.  First compute the calendar date if desired (abs(ioptn)=4).
      if (abs(ioptn) .eq. 4) then
         jmonth = (10*(ndays+31))/306 + 3
c               = offset month counter.  jmonth=4 for March, =13 for
c                 December, =14 for January, =15 for February.
         iday   = (ndays+123) - int(30.6001*jmonth)
c               = day of the month, starting with 1 on the first day
c                 of the month.
c
c        Now adjust for the fact that the year actually begins
c        on 1 January.
         if (jmonth .le. 13) then
            month = jmonth - 1
         else
            month = jmonth - 13
            iyear = iyear + 1
         endif
c
c This code handles abs(ioptn)=5, finding the day number during the year.
      else
c        ioptn=5 always returns month=1, which we set now.
         month = 1
c
c        We need to determine whether this is a leap year.
         leap = 0
         if ((jyear/4)*4 .eq. jyear) then
            leap = 1
         endif
         if ((ioptn .gt. 0)               .and.
     &       ((jyear/100)*100 .eq. jyear) .and.
     &       ((jyear/400)*400 .ne. jyear)      ) then
               leap = 0
         endif
c
c        Now find the day number "iday".
c        ndays is the number of days since the most recent 1 March,
c        so ndays=0 on 1 March.
         if (ndays .le.305) then
            iday  = ndays + 60 + leap
         else
            iday  = ndays - 305
            iyear = iyear + 1
         endif
      endif
c
c     Adjust the year if it is .le. 0, and hence BC (Before Christ).
      if (iyear .le. 0) then
         iyear = iyear - 1
      endif
c
c End the code for the last option, ioptn.
      endif
c
      return
      end


      integer function idaywk(jdayno)
c
c IDAYWK = compute the DAY of the WeeK given the Julian Day number,
c          version 1.0.
c
c Input variable
      integer  jdayno
c jdayno = Julian Day number starting at noon of the day in question.
c
c Output variable:
c idaywk = day of the week, where 0=Sunday, 1=Monday, ..., 6=Saturday.
c
c----------
c Compute the day of the week given the Julian Day number.
c You can find the Julian Day number given (day,month,year)
c using subroutine calndr.f.
c Example: For the first day of the Gregorian calendar,
c 15 October 1582, compute the Julian day number (option 3 of
c subroutine calndr) and compute the day of the week.
c     call calndr (3, 15, 10, 1582, jdayno) 
c     write(*,*) jdayno, idaywk(jdayno)
c The numbers printed should be 2299161 and 5,
c where 6 refers to Friday.
c
c Copyright (C) 1999 Jon Ahlquist.
c Issued under the second GNU General Public License.
c See www.gnu.org for details.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c If you find any errors, please notify:
c Jon Ahlquist <ahlquist@met.fsu.edu>   
c Dept of Meteorology
c Florida State University
c Tallahassee, FL 32306-4520
c 15 March 1999.
c
c-----
c Declare internal variable.
c jdSun is the Julian Day number starting at noon on any Sunday.
c I arbitrarily chose the first Sunday after Julian Day 1,
c which is Julian Day 6.
      integer  jdSun
      data     jdSun /6/
      idaywk = mod(jdayno-jdSun,7)
c If jdayno-jdSun < 0, then we are taking the modulus of a negative
c number. Fortran's built-in mod function returns a negative value
c when the argument is negative.  In that case, we adjust the result
c to a positive value.
      if (idaywk .lt. 0) idaywk = idaywk + 7
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

