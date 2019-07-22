c preprocess_code.f

c Perform a variety of preprocessing steps, like read in topography
c   and vegetation arrays, open input and output files, etc.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE PREPROCESS_CODE(topoveg_fname,const_veg_flag,
     &  vegtype,veg_z0,vegsnowdepth,fetch,xmu,C_z,h_const,
     &  wind_min,Up_const,dz_susp,ztop_susp,fall_vel,Ur_const,
     &  ro_water,ro_air,gravity,vonKarman,pi,twopio360,snow_z0,
     &  nx,ny,sum_sprec,sum_qsubl,sum_trans,sum_unload,topo,
     &  topo_land,snow_d,topoflag,snow_d_init,snow_d_init_const,
     &  soft_snow_d,met_input_fname,igrads_metfile,deltax,deltay,
     &  snowtran_output_fname,micromet_output_fname,
     &  enbal_output_fname,snowpack_output_fname,print_micromet,
     &  print_enbal,print_snowpack,print_snowtran,run_micromet,
     &  run_enbal,run_snowpack,run_snowtran,ro_snow_grid,swe_depth,
     &  sum_runoff,sum_prec,ro_snow,twolayer_flag,sum_Qcs,
     &  canopy_int,ascii_topoveg,topo_ascii_fname,icorr_factor_loop,
     &  veg_ascii_fname,undef,isingle_stn_flag,max_iter,
     &  i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag,sum_glacmelt,
     &  snow_depth,sum_d_canopy_int,corr_factor,icorr_factor_index,
     &  sum_sfcsublim,barnes_lg_domain,n_stns_used,k_stn,xmn,ymn,
     &  ro_soft_snow_old,sum_swemelt,xlat,lat_solar_flag,xlat_grid,
     &  xlon_grid,UTC_flag,dt,swe_depth_old,canopy_int_old,
     &  vegsnowd_xy,iveg_ht_flag,ihrestart_flag,i_dataassim_loop,
     &  multilayer_snowpack,max_layers,multilayer_output_fname,
     &  print_multilayer,JJ,tslsnowfall,tsls_threshold,
     &  irun_data_assim,izero_snow_date,iclear_mn,iclear_dy,
     &  xclear_hr,dy_snow,swe_lyr,ro_layer,T_old,gamma,icond_flag,
     &  curve_lg_scale_flag,curve_wt_lg,check_met_data,seaice_run,
     &  snowmodel_line_flag,xg_line,yg_line,print_user,albedo_flag)

      implicit none

      include 'snowmodel.inc'

      integer i,j,k,nx,ny,igrads_metfile,n_recs_out,iheader,
     &  isingle_stn_flag,max_iter,i_tair_flag,i_rh_flag,i_wind_flag,
     &  i_prec_flag,iter,iobs_num,n_stns_used,nveg,iveg_ht_flag,
     &  lat_solar_flag,ihrestart_flag,nstns_orig,i_dataassim_loop,
     &  multilayer_snowpack,max_layers,irun_data_assim,
     &  izero_snow_date,iclear_mn,iclear_dy,icond_flag,albedo_flag

      real ro_water,ro_air,gravity,vonKarman,snow_z0,
     &  fetch,xmu,C_z,h_const,wind_min,Up_const,check_met_data,
     &  dz_susp,ztop_susp,fall_vel,Ur_const,pi,twopio360,topoflag,
     &  snow_d_init_const,const_veg_flag,ro_snow,twolayer_flag,
     &  ascii_topoveg,undef,barnes_lg_domain,xlat,UTC_flag,dt,
     &  print_multilayer,xclear_hr,curve_lg_scale_flag,seaice_run,
     &  snowmodel_line_flag,print_user

      real topo_land(nx_max,ny_max)
      real topo(nx_max,ny_max)
      real vegtype(nx_max,ny_max)
      real xlat_grid(nx_max,ny_max)
      real xlon_grid(nx_max,ny_max)

      real snow_d(nx_max,ny_max)
      real snow_depth(nx_max,ny_max)
      real snow_d_init(nx_max,ny_max)
      real canopy_int(nx_max,ny_max)
      real swe_depth_old(nx_max,ny_max)
      real canopy_int_old(nx_max,ny_max)

      real sum_sprec(nx_max,ny_max)
      real sum_qsubl(nx_max,ny_max)
      real sum_trans(nx_max,ny_max)
      real sum_unload(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)
      real ro_soft_snow_old(nx_max,ny_max)
      real ro_snow_grid(nx_max,ny_max)
      real swe_depth(nx_max,ny_max)
      real sum_prec(nx_max,ny_max)
      real sum_runoff(nx_max,ny_max)
      real sum_Qcs(nx_max,ny_max)
      real sum_glacmelt(nx_max,ny_max)
      real sum_swemelt(nx_max,ny_max)
      real sum_d_canopy_int(nx_max,ny_max)
      real sum_sfcsublim(nx_max,ny_max)

      real vegsnowdepth(nvegtypes)
      real veg_z0(nx_max,ny_max)
      real vegsnowd_xy(nx_max,ny_max)

      real curve_wt_lg(nx_max,ny_max)

      real corr_factor(nx_max,ny_max,max_obs_dates+1)
      integer icorr_factor_index(max_time_steps)
      integer icorr_factor_loop

      integer k_stn(nx_max,ny_max,5)
      double precision xmn,ymn
      real deltax,deltay
      integer icount,iii,jjj
      double precision xg_line(nx_max,ny_max),yg_line(nx_max,ny_max)

      real run_micromet,run_enbal,run_snowpack,run_snowtran
      real print_micromet,print_enbal,print_snowpack,print_snowtran

      character*80 topoveg_fname,met_input_fname,topo_ascii_fname,
     &  veg_ascii_fname
      character*80 snowtran_output_fname,micromet_output_fname,
     &  enbal_output_fname,snowpack_output_fname,
     &  multilayer_output_fname

      integer JJ(nx_max,ny_max)
      real tslsnowfall(nx_max,ny_max)
      real tsls_threshold
      real dy_snow(nx_max,ny_max,nz_max)
      real swe_lyr(nx_max,ny_max,nz_max)
      real ro_layer(nx_max,ny_max,nz_max)
      real T_old(nx_max,ny_max,nz_max)
      real gamma(nx_max,ny_max,nz_max)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Check to see whether the maximum array dimensions defined in
c   snowmodel.inc are adequate for the simulation domain defined
c   in snowmodel.par.
      if (snowmodel_line_flag.eq.0.0) then
        if (nx+1.gt.nx_max .or. ny+1.gt.ny_max) then
          print *, 'Must increase the value of nx_max or ny_max'
          print *, '  in snowmodel.inc to be greater than nx+1'
          print *, '  and/or ny+1.'
          print *, 'nx_max = ',nx_max,'  ny_max = ',ny_max
          print *, 'nx = ',nx,'  ny = ',ny
          stop
        endif
      else
        if (nx.ge.nx_max .or. ny.ne.1 .or. ny_max.ne.2) then
          print *, 'For snowmodel_line_flag = 1.0, we suggest setting'
          print *, 'nx = number of grid cells, ny = 1, nx_max = nx+1,'
          print *, 'and ny_max = ny+1 = 2 in snowmodel.inc.'
          print *, '  The current values are:'
          print *, '    nx_max = ',nx_max,'  ny_max = ',ny_max
          print *, '    nx = ',nx,'  ny = ',ny
          stop
        endif
      endif

      if (multilayer_snowpack.eq.1) then
        if (max_layers+1.gt.nz_max) then
          print *, 'nz_max in snowmodel.inc must be at least 1 greater'
          print *, '  than max_layers in the snowmodel.inc file.  So,'
          print *, '  if you want to run the multi-layer snowpack model'
          print *, '  with a single snow layer, set nz_max=2.  If you'
          print *, '  want to run the original single-layer snowpack'
          print *, '  model, you can set nz_max=1 in snowmodel.inc.'
          print *, 'nz_max = ',nz_max
          print *, 'max_layers = ',max_layers
          stop
        endif
      endif

      if (max_iter.gt.max_time_steps) then
        print *, 'Must increase the value of max_time_steps'
        print *, '  in snowmodel.inc to be greater than max_iter.'
        print *, 'max_time_steps = ',max_time_steps
        print *, 'max_iter = ',max_iter
        stop
      endif

c If running the concatenated configuration of the model, check to
c   make sure the rest of the model is configured correctly.
      if (snowmodel_line_flag.eq.1.0) then
        if (run_snowtran.eq.1.0) then
          print *, 'You cannot run snowmodel_line_flag = 1.0 with'
          print *, 'run_snowtran = 1.0'
        stop
        endif
        if (barnes_lg_domain.eq.0.0) then
          print *, 'If snowmodel_line_flag = 1.0, then you must run'
          print *, 'the model with barnes_lg_domain = 1.0'
        stop
        endif
      endif

c Make sure the time since last snowfall treshold is not less
c   than the model time step.
      if (multilayer_snowpack.eq.1) then
        if (tsls_threshold.lt.dt/3600.0) then
          print *,'Need to correct tsls_threshold to'
          print *,'  be >= dt (in hours).'
          stop
        endif
      endif

c Check the model dt value, and send an error message if dt < 3600.
      if (dt.lt.3600.0) then
        print *, 'You must modify the hour fraction calculation'
        print *, '  in get_model_time subroutine to handle'
        print *, '  dt values less that 3600.0 seconds.'
        print *, 'dt = ',dt
        stop
      endif

c Define the date on which the snow arrays will be zeroed out.
      iclear_mn = izero_snow_date / 10000
      iclear_dy = (izero_snow_date - iclear_mn * 10000) / 100
      xclear_hr =
     &  real((izero_snow_date - iclear_mn * 10000) - iclear_dy * 100)

c Check to see whether there is enough snow layers to calculate
c   conductive surface fluxes.
      if (icond_flag.eq.1) then
        if (multilayer_snowpack.eq.0 .or. max_layers.lt.2) then
          print *,'If icond_flag = 1, then multilayer_snowpack = 1'
          print *,'  and max_layers >= 2.'
          stop
        endif
      endif

c Read in the topography array.
      if (ascii_topoveg.eq.0.0) then

        open (unit=37,file=topoveg_fname,
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (37,rec=1) ((topo_land(i,j),i=1,nx),j=1,ny)

      elseif (ascii_topoveg.eq.1.0) then

c Read off the header lines.  I will assume that all of this
c   information was input in the .par file correctly.
        open (37,file=topo_ascii_fname,form='formatted')
        iheader = 6
        do k=1,iheader
          read (37,*)
        enddo
c Read the data in as real numbers, and do the yrev.
        do j=ny,1,-1
          read (37,*) (topo_land(i,j),i=1,nx)
        enddo

      endif

c If vegetation data is not available on the topography grid,
c   define the vegetation to be constant.
      if (const_veg_flag.ne.0.0) then
        do i=1,nx
          do j=1,ny
            vegtype(i,j) = const_veg_flag
          enddo
        enddo

c Read in the vegetation array.
      else

        if (ascii_topoveg.eq.0.0) then
          read (37,rec=2) ((vegtype(i,j),i=1,nx),j=1,ny)

        elseif (ascii_topoveg.eq.1.0) then

c Read off the header lines.  I will assume that all of this
c   information was input in the .par file correctly.
          open (38,file=veg_ascii_fname,form='formatted')
          do k=1,iheader
            read (38,*)
          enddo
c Read the data in as real numbers, and do the yrev.
          do j=ny,1,-1
            read (38,*) (vegtype(i,j),i=1,nx)
          enddo

        endif

      endif

c Now that we have read in the topo and veg data arrays, check
c   whether all of the values look like valid numbers.
      do i=1,nx
        do j=1,ny
          if (vegtype(i,j).lt.1.0 .or. vegtype(i,j).gt.30.0) then
            print *, 'Found Invalid Vegetation-Type Value'
            print *, '     Value =',vegtype(i,j),'  at ',i,j
            stop
          endif

          if (topo_land(i,j).lt.0.0 .or. topo_land(i,j).gt.9000.0) then
            print *, 'Found Invalid Topography Value'
            print *, '     Value =',topo_land(i,j),'  at ',i,j
            stop
          endif
        enddo
      enddo

c J.PFLUG
c provide the opportunity to read in snow and rain precipitation data
      if (i_prec_flag.eq.-1.0) then
        open (99,file='met/TUO_wind_fp.txt',
     &    form='formatted')
        rewind(99)
      endif

      if (albedo_flag.eq.1.0) then
        open (98,file='extra_met/albedo.dat',
     &    form='formatted')
        rewind (98)
      endif
c END J.PFLUG

c Fill the the vegetation snow-holding depth array for vegetation
c   types 1 through 24 (types 25 through 30 were filled from the
c   .par file.
      call fill_veg_shd(nvegtypes,vegsnowdepth)

c Use vegsnowdepth to fill the 2D spatial array, or read in the
c   user-provided file of the vegetation heights (in cm).
      if (iveg_ht_flag.eq.-1) then

        open (191,file='topo_veg/veg_ht.gdat',
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (191,rec=1) ((vegsnowd_xy(i,j),i=1,nx),j=1,ny)
        close(191)

c Convert from cm to m.
        do i=1,nx
          do j=1,ny
            vegsnowd_xy(i,j) = vegsnowd_xy(i,j) / 100.0
          enddo
        enddo

      elseif (iveg_ht_flag.eq.1) then

        open (191,file='topo_veg/veg_ht.asc',form='formatted')
        iheader = 6
        do k=1,iheader
          read (191,*)
        enddo
c Read the data in as real numbers, and do the yrev.
        do j=ny,1,-1
          read (191,*) (vegsnowd_xy(i,j),i=1,nx)
        enddo
        close(191)

c Convert from cm to m.
        do i=1,nx
          do j=1,ny
            vegsnowd_xy(i,j) = vegsnowd_xy(i,j) / 100.0
          enddo
        enddo

      elseif (iveg_ht_flag.eq.0) then

        do i=1,nx
          do j=1,ny
            nveg = nint(vegtype(i,j))
            vegsnowd_xy(i,j) = vegsnowdepth(nveg)
          enddo
        enddo

      endif

c Define the roughness lengths for each of the vegetation types.
c   Note that this is pretty arbitrary, because these values are
c   really only used when there is no blowing snow, and thus have
c   no impact on the simulation except to provide a non-zero value
c   for any such parts of the domain.
      do i=1,nx
        do j=1,ny
          veg_z0(i,j) = 0.25 * vegsnowd_xy(i,j)
        enddo
      enddo

c Read in the large-scale curvature weighting array, if the run
c   requires it.
      if (curve_lg_scale_flag.eq.1.0) then
        open (444,file='extra_met/large_curvature_wt.gdat',
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (444,rec=1) ((curve_wt_lg(i,j),i=1,nx),j=1,ny)
        close (444)
      endif

c If this is a sea ice run, open the sea ice concentration file.
      if (seaice_run.ne.0.0) then
        open (445,file='seaice/ice_conc.gdat',
     &    form='unformatted',access='direct',recl=4*nx*ny)
      endif

c Check to make sure that if you are running SnowTran-3D and the
c   EnBal and SnowPack models, you have also set the flag to run
c   the SnowTran-3D two-layer submodel.
      if (run_enbal.eq.1.0 .and. run_snowpack.eq.1.0 .and.
     &  run_snowtran.eq.1.0) then
        if (twolayer_flag.ne.1.0) then
          print *, 'For SnowTran-3D with EnBal and SnowPack,'
          print *, '  twolayer_flag must equal 1.0'
          stop
        endif
      endif

c Check to see that the defined domain is large enough to be
c   running SnowTran-3D.
      if (run_snowtran.eq.1.0) then
        if (nx.lt.3 .or. ny.lt.3) then
          print *, 'To run SnowTran-3D, nx and ny must both be 3'
          print *, '  or greater (see SnowTran-3D code/notes)'
          stop
        endif
      endif

c Check to see whether the model is configured correctly to be
c   running the multi-layer snow model.
      if (multilayer_snowpack.eq.1) then
        if (irun_data_assim.eq.1 .or. ihrestart_flag.ne.-2 .or.
     &    snow_d_init_const.ne.0.0) then
          print *, 'The multi-layer snowpack model requires:'
          print *, '  irun_data_assim = 0'
          print *, '  ihrestart_flag = -2'
          print *, '  snow_d_init_const = 0.0'
c          stop
        endif
      endif

c Get a collection of constants that are not generally modified.
      call constants(fetch,xmu,C_z,h_const,wind_min,Up_const,
     &  dz_susp,ztop_susp,fall_vel,Ur_const,ro_water,ro_air,
     &  gravity,vonKarman,pi,twopio360,snow_z0)

c Run a check to see if SnowTran-3D is being run with grid
c   increments that are too large.
      if (deltax.gt.500.0 .or. deltay.gt.500.0) then
        if (run_snowtran.eq.1.0) then
          print *
          print *, '!!! deltax,y should not be greater than 500 m'
          print *, '!!!    if you are also running SnowTran-3D'
          print *
          stop
        endif
      endif

c Initialize the summing arrays, and define the initial snow-depth
c   distributions.
      call initialize(nx,ny,sum_sprec,sum_qsubl,sum_trans,
     &  sum_unload,topo,topo_land,snow_d,topoflag,snow_d_init,
     &  snow_d_init_const,soft_snow_d,ro_water,sum_sfcsublim,
     &  ro_snow_grid,swe_depth,sum_runoff,sum_prec,ro_snow,
     &  sum_Qcs,canopy_int,sum_glacmelt,snow_depth,sum_d_canopy_int,
     &  ro_soft_snow_old,sum_swemelt,swe_depth_old,canopy_int_old,
     &  ihrestart_flag,i_dataassim_loop,max_iter,corr_factor,
     &  icorr_factor_index,JJ,tslsnowfall,tsls_threshold,dy_snow,
     &  swe_lyr,ro_layer,T_old,gamma)

c Check to see whether the data assimilation has been configured
c   correctly.
      if (irun_data_assim.eq.1) then

c Check to see whether the required output files will be created.
        if (print_user.ne.1.0) then
          print *, 'For a data assimilation run print_user must = 1.0'
          stop
        endif

c Check to see whether the corr_factor array is defined to be large
c   enough to do the assimilation.  It looks like I can't do that
c   here because I don't know the name of the swe-obs input file yet.
c   Therefore I will do a check when dataassim_user.f is called.  So
c   in the mean time I will just print a message to the screen telling
c   the user to check this.
        print *
        print *, 'For a DATA ASSIMILATION RUN, MAX_OBS_DATES must be'
        print *, 'defined in SNOWMODEL.INC to be greater than the'
        print *, 'number of observation dates in the entire simulation'
        print *, '+ (plus) the number of years in the simulation.  For'
        print *, 'example, for a 6-year simulation with two observation'
        print *, 'dates in each year, you would set max_obs_dates to be'
        print *, 'at least = 18.'
        print *

      endif

c Initialize the precipitation factor for the first iteration to
c   equal 1.0.
      if (icorr_factor_loop.eq.1) then
        do iobs_num=1,max_obs_dates+1
          do j=1,ny
            do i=1,nx
              corr_factor(i,j,iobs_num) = 1.0
            enddo
          enddo
        enddo
        do iter=1,max_iter
          icorr_factor_index(iter) = 1
        enddo
      endif

c Read or build the latitude array that will be used to do the
c   latitude weighting when calculating incoming solar radiation.
      if (lat_solar_flag.eq.-1) then

        open (91,file='extra_met/grid_lat.gdat',
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (91,rec=1) ((xlat_grid(i,j),i=1,nx),j=1,ny)
        close(91)

      elseif (lat_solar_flag.eq.1) then

        open (91,file='extra_met/grid_lat.asc',form='formatted')
        iheader = 6
        do k=1,iheader
          read (91,*)
        enddo
c Read the data in as real numbers, and do the yrev.
        do j=ny,1,-1
          read (91,*) (xlat_grid(i,j),i=1,nx)
        enddo
        close(91)

      elseif (lat_solar_flag.eq.0) then

c Print an error if the y-domain is big enough to have important
c   solar radiation differences from south to north.
        if (ny*deltay.gt.500000.0) then
          print *
          print *,'YOUR DOMAIN IS PRETTY BIG TO NOT ACCOUNT FOR'
          print *,'  SOLAR RADIATION VARIATIONS WITH LATITUDE'
          print *,' see the "lat_solar_flag" in snowmodel.par'
          print *
          stop
        endif

        do i=1,nx
          do j=1,ny
            xlat_grid(i,j) = xlat
          enddo
        enddo

      endif

c Read or build the lonitude array that will be used to do the
c   lonitude influence when calculating incoming solar radiation.
      if (UTC_flag.eq.-1.0) then

        open (91,file='extra_met/grid_lon.gdat',
     &    form='unformatted',access='direct',recl=4*nx*ny)
        read (91,rec=1) ((xlon_grid(i,j),i=1,nx),j=1,ny)
        close(91)

      elseif (UTC_flag.eq.1.0) then

        open (91,file='extra_met/grid_lon.asc',form='formatted')
        iheader = 6
        do k=1,iheader
          read (91,*)
        enddo
c Read the data in as real numbers, and do the yrev.
        do j=ny,1,-1
          read (91,*) (xlon_grid(i,j),i=1,nx)
        enddo
        close(91)

      elseif (UTC_flag.eq.0.0) then

c Print an error if the x-domain is big enough to have important
c   solar radiation differences from east to west.
        if (nx*deltax.gt.500000.0) then
          print *
          print *,'YOUR DOMAIN IS PRETTY BIG TO NOT ACCOUNT FOR'
          print *,'  SOLAR RADIATION VARIATIONS WITH LONITUDE'
          print *,'    see the "UTC_flag" in snowmodel.par'
          print *
          stop
        endif

      endif

c Open the MicroMet station data input file.
      if (igrads_metfile.eq.1) then
        open(20,file=met_input_fname,form='unformatted',
     &    access='direct',recl=4*13)
      else
        open (20,file=met_input_fname,form='formatted')
      endif

c Run a check to see whether there are any time slices with no
c   valid data.
      if (check_met_data.eq.1.0) then
        print *
        print *,'Checking for sufficient met forcing data to'
        print *,'  complete the model simulation.  This may'
        print *,'  take a while, depending on how big your met'
        print *,'  input file is.'
        print *
        call met_data_check(undef,isingle_stn_flag,igrads_metfile,
     &    max_iter,i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag)
      endif

c If the concatenated configuration of the model is used, read
c   in the x and y coordinates for the concatenated grid cells.
      if (snowmodel_line_flag.eq.1.0) then
        open (1331,file='extra_met/snowmodel_line_pts.dat')
        do j=1,ny
          do i=1,nx
            read (1331,*) icount,iii,jjj,xg_line(i,j),yg_line(i,j)
          enddo
        enddo
        close (1331)
      endif

c If the large-domain barnes oi scheme is used, generate the
c   nearest-station indexing array.
      if (barnes_lg_domain.eq.1.0) then
        print *
        print *,'You are running the large-domain Barnes oi scheme'
        print *,'  This requires:'
        print *,'  1) no missing data for the fields of interest'
        print *,'  2) no missing stations during the simulation' 
        print *,'  3) met file must list stations in the same order'
        print *,'  4) the number of nearest stations used is 5 or less'
        print *,'  5)  **** no error checking for this is done ****'
        print *
        print *,'Generating nearest-station index.  Be patient.'
        print *
        if (n_stns_used.gt.5 .or. n_stns_used.lt.1) then
          print *,'invalid n_stns_used value'
          stop
        endif
        call get_nearest_stns_1(nx,ny,xmn,ymn,deltax,deltay,
     &    n_stns_used,k_stn,snowmodel_line_flag,xg_line,yg_line)
      endif

c If this is a history restart run, advance the micromet input
c   file to the restart time.
      if (ihrestart_flag.ge.0) then
        if (igrads_metfile.eq.0) then
          do iter=1,ihrestart_flag
            if (isingle_stn_flag.eq.1) then
              nstns_orig = 1
            else
              read(20,*) nstns_orig
            endif
            do k=1,nstns_orig
              read(20,*)
            enddo
          enddo
        endif
      endif

c Open the files to be used to store model output.
c   For MicroMet.
      if (run_micromet.eq.1.0 .and. print_micromet.eq.1.0) then
        n_recs_out = 8
        open (81,file=micromet_output_fname,
     &    form='unformatted',access='direct',recl=4*n_recs_out*nx*ny)
      endif

c   For EnBal.
      if (run_enbal.eq.1.0 .and. print_enbal.eq.1.0) then
        n_recs_out = 10
        open (82,file=enbal_output_fname,
     &    form='unformatted',access='direct',recl=4*n_recs_out*nx*ny)
      endif

c   For SnowPack.
      if (run_snowpack.eq.1.0 .and. print_snowpack.eq.1.0) then
        n_recs_out = 16
        open (83,file=snowpack_output_fname,
     &    form='unformatted',access='direct',recl=4*n_recs_out*nx*ny)
      endif

c   For SnowTran-3D.
      if (run_snowtran.eq.1.0 .and. print_snowtran.eq.1.0) then
        n_recs_out = 7
        open (84,file=snowtran_output_fname,
     &    form='unformatted',access='direct',recl=4*n_recs_out*nx*ny)
      endif

c   For Multi-Layer SnowPack.
      if (run_snowpack.eq.1.0 .and. multilayer_snowpack.eq.1 .and.
     &  print_multilayer.eq.1.0) then
        open (85,file=multilayer_output_fname,
     &    form='unformatted',access='direct',
     &    recl=4*(4*nx*ny+6*nx*ny*nz_max))
      endif

c    For random precipitation error
      if (irun_data_assim.eq.1.and.icorr_factor_loop.eq.1) then
        open(11,file='data/randomError.txt')
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initialize(nx,ny,sum_sprec,sum_qsubl,sum_trans,
     &  sum_unload,topo,topo_land,snow_d,topoflag,snow_d_init,
     &  snow_d_init_const,soft_snow_d,ro_water,sum_sfcsublim,
     &  ro_snow_grid,swe_depth,sum_runoff,sum_prec,ro_snow,
     &  sum_Qcs,canopy_int,sum_glacmelt,snow_depth,sum_d_canopy_int,
     &  ro_soft_snow_old,sum_swemelt,swe_depth_old,canopy_int_old,
     &  ihrestart_flag,i_dataassim_loop,max_iter,corr_factor,
     &  icorr_factor_index,JJ,tslsnowfall,tsls_threshold,dy_snow,
     &  swe_lyr,ro_layer,T_old,gamma)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,ihrestart_flag,i_dataassim_loop,max_iter,k

      real topoflag,snow_d_init_const,ro_water,ro_snow
      real sum_sprec(nx_max,ny_max)
      real sum_qsubl(nx_max,ny_max)
      real sum_trans(nx_max,ny_max)
      real sum_unload(nx_max,ny_max)
      real topo(nx_max,ny_max)
      real topo_land(nx_max,ny_max)
      real snow_d(nx_max,ny_max)
      real snow_depth(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)
      real ro_soft_snow_old(nx_max,ny_max)
      real ro_snow_grid(nx_max,ny_max)
      real swe_depth(nx_max,ny_max)
      real sum_prec(nx_max,ny_max)
      real sum_runoff(nx_max,ny_max)
      real sum_Qcs(nx_max,ny_max)
      real canopy_int(nx_max,ny_max)
      real sum_glacmelt(nx_max,ny_max)
      real sum_swemelt(nx_max,ny_max)
      real sum_d_canopy_int(nx_max,ny_max)
      real sum_sfcsublim(nx_max,ny_max)
      real snow_d_init(nx_max,ny_max)
      real swe_depth_old(nx_max,ny_max)
      real canopy_int_old(nx_max,ny_max)

      integer JJ(nx_max,ny_max)
      real tslsnowfall(nx_max,ny_max)
      real tsls_threshold
      real dy_snow(nx_max,ny_max,nz_max)
      real swe_lyr(nx_max,ny_max,nz_max)
      real ro_layer(nx_max,ny_max,nz_max)
      real T_old(nx_max,ny_max,nz_max)
      real gamma(nx_max,ny_max,nz_max)

      integer icorr_factor_index(max_time_steps)
      real corr_factor(nx_max,ny_max,max_obs_dates+1)

      if (ihrestart_flag.ge.0) then

c Read in the saved data.
        CALL HRESTART_READ(nx,ny,snow_d,snow_depth,
     &    canopy_int,soft_snow_d,ro_snow_grid,swe_depth,
     &    ro_soft_snow_old,snow_d_init,swe_depth_old,
     &    canopy_int_old,topo,sum_sprec,ihrestart_flag,
     &    i_dataassim_loop)

        if (i_dataassim_loop.lt.0.0) then
          CALL HRESTART_READ_DA(nx,ny,max_iter,corr_factor,
     &      icorr_factor_index,i_dataassim_loop)
        endif

        do i=1,nx
          do j=1,ny
c Fill the summing arrays.
            sum_runoff(i,j) = 0.0
            sum_prec(i,j) = 0.0
c           sum_sprec(i,j) = 0.0
            sum_qsubl(i,j) = 0.0
            sum_trans(i,j) = 0.0
            sum_unload(i,j) = 0.0
            sum_Qcs(i,j) = 0.0
            sum_glacmelt(i,j) = 0.0
            sum_swemelt(i,j) = 0.0
            sum_d_canopy_int(i,j) = 0.0
            sum_sfcsublim(i,j) = 0.0

c Define the initial snow-depth distributions.
c           snow_d_init(i,j) = snow_d_init_const
c           snow_d(i,j) = snow_d_init(i,j)
c           snow_depth(i,j) = snow_d_init(i,j)
c           canopy_int(i,j) = 0.0
c           soft_snow_d(i,j) = snow_d(i,j)
c           ro_snow_grid(i,j) = ro_snow
c           swe_depth(i,j) = snow_d(i,j) * ro_snow_grid(i,j) / ro_water
c           ro_soft_snow_old(i,j) = 50.0
c           swe_depth_old(i,j) = swe_depth(i,j)
c           canopy_int_old(i,j) = canopy_int(i,j)
          enddo
        enddo

      else

        do i=1,nx
          do j=1,ny
c Fill the summing arrays.
            sum_runoff(i,j) = 0.0
            sum_prec(i,j) = 0.0
            sum_sprec(i,j) = 0.0
            sum_qsubl(i,j) = 0.0
            sum_trans(i,j) = 0.0
            sum_unload(i,j) = 0.0
            sum_Qcs(i,j) = 0.0
            sum_glacmelt(i,j) = 0.0
            sum_swemelt(i,j) = 0.0
            sum_d_canopy_int(i,j) = 0.0
            sum_sfcsublim(i,j) = 0.0

c Define the initial snow-depth distributions.
            snow_d_init(i,j) = snow_d_init_const
            snow_d(i,j) = snow_d_init(i,j)
            snow_depth(i,j) = snow_d_init(i,j)
            canopy_int(i,j) = 0.0
            soft_snow_d(i,j) = snow_d(i,j)
            ro_snow_grid(i,j) = ro_snow
            swe_depth(i,j) = snow_d(i,j) * ro_snow_grid(i,j) / ro_water
            ro_soft_snow_old(i,j) = 50.0
            swe_depth_old(i,j) = swe_depth(i,j)
            canopy_int_old(i,j) = canopy_int(i,j)

c Initialize the multi-layer snowpack arrays.
            JJ(i,j) = 0
            tslsnowfall(i,j) = tsls_threshold
          enddo
        enddo

        do i=1,nx
          do j=1,ny
            do k=1,nz_max
              dy_snow(i,j,k) = 0.0
              swe_lyr(i,j,k) = 0.0
              ro_layer(i,j,k) = ro_snow
              T_old(i,j,k) = 273.16
              gamma(i,j,k) = 0.138 - 1.01 * (ro_layer(i,j,k)/1000.0) +
     &          3.233 * (ro_layer(i,j,k)/1000.0)**2
            enddo
          enddo
        enddo

        if (topoflag.eq.1.0) then
          do i=1,nx
            do j=1,ny
              topo(i,j) = topo_land(i,j) + snow_d(i,j)
            enddo
          enddo
        elseif (topoflag.eq.0.0) then
          do i=1,nx
            do j=1,ny
              topo(i,j) = topo_land(i,j)
            enddo
          enddo
        endif

      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine constants(fetch,xmu,C_z,h_const,wind_min,Up_const,
     &  dz_susp,ztop_susp,fall_vel,Ur_const,ro_water,ro_air,
     &  gravity,vonKarman,pi,twopio360,snow_z0)

      implicit none

      real fetch,xmu,C_z,h_const,wind_min,Up_const,
     &  dz_susp,ztop_susp,fall_vel,Ur_const,ro_water,ro_air,
     &  gravity,vonKarman,pi,twopio360,snow_z0

c These constants are not generally modified for a particular model
c   run.

c Snow surface roughness length.
      snow_z0 = 0.001

c Constants related to surface shear stress and saltation
c   transport.
      fetch = 500.0
      xmu = 3.0
      C_z = 0.12
      h_const = 1.6
      wind_min = 4.0

c Constants related to suspended snow profile.
      Up_const = 2.8
      dz_susp = 0.20
      ztop_susp = 2.0
      fall_vel = 0.3
      Ur_const = 0.5

c General constants.
      ro_water = 1000.0
      ro_air = 1.275
      gravity = 9.81
      vonKarman = 0.4
      pi = 2.0 * acos(0.0)
      twopio360 = 2.0 * pi / 360.0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fill_veg_shd(nvegtypes,vegsnowdepth)

      implicit none

      integer k,nvegtypes,nvegtypes_fixed

      parameter (nvegtypes_fixed=24)

      real vegsnowdepth(nvegtypes),vegsnowdepth_fixed(nvegtypes_fixed)

c Fill the the vegetation snow-holding depth array for
c   vegetation types 1 through 24 (types 25 through 30 were filled
c   from the .par file.
c
c The following summary was taken from the .par file.
c
c The vegetation types are assumed to range from 1 through 30.  The
c   last 6 types are available to be user-defined vegetation types
c   and vegetation snow-holding depth.  The first 24 vegetation
c   types, and the associated vegetation snow-holding depth
c   (meters), are hard-coded to be:
c
c code description           veg_shd  example                    class
c
c  1  coniferous forest       15.00  spruce-fir/taiga/lodgepole  forest
c  2  deciduous forest        12.00  aspen forest                forest
c  3  mixed forest            14.00  aspen/spruce-fir/low taiga  forest
c  4  scattered short-conifer  8.00  pinyon-juniper              forest
c  5  clearcut conifer         4.00  stumps and regenerating     forest
c 
c  6  mesic upland shrub       0.50  deeper soils, less rocky    shrub
c  7  xeric upland shrub       0.25  rocky, windblown soils      shrub
c  8  playa shrubland          1.00  greasewood, saltbush        shrub
c  9  shrub wetland/riparian   1.75  willow along streams        shrub
c 10  erect shrub tundra       0.65  arctic shrubland            shrub
c 11  low shrub tundra         0.30  low to medium arctic shrubs shrub
c 
c 12  grassland rangeland      0.15  graminoids and forbs        grass
c 13  subalpine meadow         0.25  meadows below treeline      grass
c 14  tundra (non-tussock)     0.15  alpine, high arctic         grass
c 15  tundra (tussock)         0.20  graminoid and dwarf shrubs  grass
c 16  prostrate shrub tundra   0.10  graminoid dominated         grass
c 17  arctic gram. wetland     0.20  grassy wetlands, wet tundra grass
c 
c 18  bare                     0.01                              bare
c
c 19  water/possibly frozen    0.01                              water
c 20  permanent snow/glacier   0.01                              water
c 
c 21  residential/urban        0.01                              human
c 22  tall crops               0.40  e.g., corn stubble          human
c 23  short crops              0.25  e.g., wheat stubble         human
c 24  ocean                    0.01                              water
c
c 25  user defined (see below)
c 26  user defined (see below)
c 27  user defined (see below)
c 28  user defined (see below)
c 29  user defined (see below)
c 30  user defined (see below)
c
c Define the vegetation snow-holding depth (meters) for each
c   of the user-defined vegetation types.  The numbers in the
c   list below correspond to the vegetation-type numbers in the
c   vegetation-type data array (veg type 25.0 -> veg_shd_25).  Note
c   that even if these are not used, they cannot be commented out
c   or you will get an error message.
c     veg_shd_25 = 0.10
c     veg_shd_26 = 0.10
c     veg_shd_27 = 0.10
c     veg_shd_28 = 0.10
c     veg_shd_29 = 0.10
c     veg_shd_30 = 0.10

      data vegsnowdepth_fixed/15.00, 12.00, 14.00,  8.00,  4.00,
     &                         0.50,  0.25,  1.00,  1.75,  0.65,  0.30,
     &                         0.15,  0.25,  0.15,  0.20,  0.10,  0.20,
     &                         0.01,  0.01,  0.01,  0.01,  0.40,  0.25,
     &                         0.01/

      do k=1,nvegtypes_fixed
        vegsnowdepth(k) = vegsnowdepth_fixed(k)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine met_data_check(undef,isingle_stn_flag,igrads_metfile,
     &  max_iter,i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag)

      implicit none

      include 'snowmodel.inc'

      integer iyr,imo,idy      ! year, month, and day of data
      real xhr                 ! decimal hour
      integer idstn            ! station id number

      integer k,nstns_orig,isingle_stn_flag,igrads_metfile,iter,
     &  n_good_Tair,n_good_rh,n_good_wspd,n_good_wdir,n_good_prec,
     &  n_notgood_vars,i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag,
     &  max_iter

      real Tair_orig(nstns_max),rh_orig(nstns_max)
      real winddir_orig(nstns_max),windspd_orig(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)
      real elev_orig(nstns_max),prec_orig(nstns_max)
      real undef               ! undefined value
      real elevation_flag

      n_notgood_vars = 0
      elevation_flag = 0.0

      do iter=1,max_iter

        n_good_Tair = 0
        n_good_rh = 0
        n_good_wspd = 0
        n_good_wdir = 0
        n_good_prec = 0

        if (igrads_metfile.eq.1) then
          nstns_orig = 1
        else
          if (isingle_stn_flag.eq.1) then
            nstns_orig = 1
          else
            read(20,*) nstns_orig
          endif
        endif

        if (nstns_orig.gt.nstns_max) then
          print *, 'The number of met stations in your MicroMet'
          print *, 'input file exceeds nstns_max in snowmodel.inc.'
          print *, 'This occurs at iter =',iter
          stop
        endif

        do k=1,nstns_orig

          if (igrads_metfile.eq.1) then
            read(20,rec=iter) iyr,imo,idy,xhr,idstn,xstn_orig(k),
     &        ystn_orig(k),elev_orig(k),Tair_orig(k),rh_orig(k),
     &        windspd_orig(k),winddir_orig(k),prec_orig(k)
          else
            read(20,*) iyr,imo,idy,xhr,idstn,xstn_orig(k),
     &        ystn_orig(k),elev_orig(k),Tair_orig(k),rh_orig(k),
     &        windspd_orig(k),winddir_orig(k),prec_orig(k)
          endif

c Count the good values at this time.
          if (Tair_orig(k).ne.undef) n_good_Tair = n_good_Tair + 1
          if (rh_orig(k).ne.undef) n_good_rh = n_good_rh + 1
          if (windspd_orig(k).ne.undef) n_good_wspd = n_good_wspd + 1
          if (winddir_orig(k).ne.undef) n_good_wdir = n_good_wdir + 1
          if (prec_orig(k).ne.undef) n_good_prec = n_good_prec + 1

          if (elev_orig(k).lt.0.0) then
            elevation_flag = 1.0
            print *,'elevation = ',elev_orig(k),'  for stn id = ',idstn
          endif

        enddo

c Check to see whether there are any variables with no valid data
c   at this time slice.
        if (n_good_Tair.eq.0 .and. i_tair_flag.eq.1) then
          n_notgood_vars = n_notgood_vars + 1
          print *,'no good Tair data at           ',iyr,imo,idy,xhr
        endif

        if (n_good_rh.eq.0 .and. i_rh_flag.eq.1) then
          n_notgood_vars = n_notgood_vars + 1
          print *,'no good rh data at             ',iyr,imo,idy,xhr
        endif

        if (n_good_wspd.eq.0 .and. i_wind_flag.eq.1) then
          n_notgood_vars = n_notgood_vars + 1
          print *,'no good wind speed data at     ',iyr,imo,idy,xhr
        endif

        if (n_good_wdir.eq.0 .and. i_wind_flag.eq.1) then
          n_notgood_vars = n_notgood_vars + 1
          print *,'no good wind direction data at ',iyr,imo,idy,xhr
        endif

        if (n_good_prec.eq.0 .and. i_prec_flag.eq.1) then
          n_notgood_vars = n_notgood_vars + 1
          print *,'no good precipitation data at  ',iyr,imo,idy,xhr
        endif

      enddo

      if (n_notgood_vars.gt.0) then
        print *
        print *,' FOUND TIMES WITH NO VALID MET OBSERVATIONS'
        print *,'NEED TO CORRECT THE PROBLEM BEFORE CONTINUING'
        stop
      endif

      if (elevation_flag.eq.1.0) then
        print *
        print *,' FOUND A NEGATIVE OR UNDEFINED STATION ELEVATION.'
        print *,'STATION ELEVATIONS CANNOT BE UNDEFINED, BUT THEY.'
        print *,'CAN BE NEGATIVE FOR A PLACE LIKE DEATH VALLEY.'
        print *,'YOU NEED TO CORRECT ANY PROBLEMS BEFORE CONTINUING.'
        stop
      endif

      if (igrads_metfile.eq.0) rewind (20)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_nearest_stns_1(nx,ny,xmn,ymn,deltax,deltay,
     &  n_stns_used,k_stn,snowmodel_line_flag,xg_line,yg_line)

      implicit none

      include 'snowmodel.inc'

      double precision xstn(nstns_max)
      double precision ystn(nstns_max)
      double precision dsq(nstns_max)
      double precision xg_line(nx_max,ny_max),yg_line(nx_max,ny_max)
      real snowmodel_line_flag

      double precision xg,yg,xmn,ymn,dist_min
      real deltax,deltay,x1,x2,x3,x4,x5,x6,x7

      integer i,j,k,kk,nstns,n_stns_used,nx,ny,i1,i2,i3,i4
      integer k_stn(nx_max,ny_max,5)

c Read the station information for the first (and all) time step(s).
      read(20,*) nstns
      do k=1,nstns
        read(20,*) i1,i2,i3,x1,i4,xstn(k),ystn(k),
     &    x2,x3,x4,x5,x6,x7
      enddo
      rewind (20)

      do j=1,ny
        do i=1,nx

c xcoords of grid nodes at index i,j
c ycoords of grid nodes at index i,j
          if (snowmodel_line_flag.eq.1.0) then
            xg = xg_line(i,j)
            yg = yg_line(i,j)
          else
            xg = xmn + deltax * (real(i) - 1.0)
            yg = ymn + deltay * (real(j) - 1.0)
          endif

c Loop through all of the stations, calculating the distance
c   between the current grid point and each of the stations.
          do k=1,nstns
            dsq(k) = (xg - xstn(k))**2 + (yg - ystn(k))**2
          enddo

c Loop through each of the station distances and find the
c   stations closest to the grid point in question.
          do kk=1,n_stns_used
            dist_min = 1.0e30
            do k=1,nstns
              if (dsq(k).le.dist_min) then
                k_stn(i,j,kk) = k
                dist_min = dsq(k)
              endif
            enddo

c Eliminate the last found minimum from the next search by making
c   its distance a big number.
            dsq(k_stn(i,j,kk)) = 1.0e30
          enddo

        enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine HRESTART_SAVE(nx,ny,iter,snow_d,snow_depth,
     &  canopy_int,soft_snow_d,ro_snow_grid,swe_depth,
     &  ro_soft_snow_old,snow_d_init,swe_depth_old,
     &  canopy_int_old,topo,sum_sprec,icorr_factor_loop,
     &  max_iter)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,iter,icorr_factor_loop,max_iter
      real snow_d(nx_max,ny_max)
      real snow_depth(nx_max,ny_max)
      real canopy_int(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)
      real ro_snow_grid(nx_max,ny_max)
      real swe_depth(nx_max,ny_max)
      real ro_soft_snow_old(nx_max,ny_max)
      real snow_d_init(nx_max,ny_max)
      real swe_depth_old(nx_max,ny_max)
      real canopy_int_old(nx_max,ny_max)
      real topo(nx_max,ny_max)
      real sum_sprec(nx_max,ny_max)

      character*18 name1
      character*1 name2
      character*5 name3
      character*5 niter
      character*1 iloop
      character*30 fname

c Build the file name so it includes the interation number.
      name1 = 'hrestart/hrestart_'
      name2 = '_'
      name3 = '.gdat'

      write(niter,'(i5.5)') iter
      write(iloop,'(i1.1)') icorr_factor_loop
      fname = name1//niter//name2//iloop//name3
c Save the data.
      open(151,file=fname,
     &  form='unformatted',access='direct',recl=4*nx*ny)

      write (151,rec=1) ((snow_d(i,j),i=1,nx),j=1,ny)
      write (151,rec=2) ((snow_depth(i,j),i=1,nx),j=1,ny)
      write (151,rec=3) ((canopy_int(i,j),i=1,nx),j=1,ny)
      write (151,rec=4) ((soft_snow_d(i,j),i=1,nx),j=1,ny)
      write (151,rec=5) ((ro_snow_grid(i,j),i=1,nx),j=1,ny)
      write (151,rec=6) ((swe_depth(i,j),i=1,nx),j=1,ny)
      write (151,rec=7) ((ro_soft_snow_old(i,j),i=1,nx),j=1,ny)
      write (151,rec=8) ((snow_d_init(i,j),i=1,nx),j=1,ny)
      write (151,rec=9) ((swe_depth_old(i,j),i=1,nx),j=1,ny)
      write (151,rec=10) ((canopy_int_old(i,j),i=1,nx),j=1,ny)
      write (151,rec=11) ((topo(i,j),i=1,nx),j=1,ny)
      write (151,rec=12) ((sum_sprec(i,j),i=1,nx),j=1,ny)

      close (151)
c Save a copy that can be used as the initial condition for the
c   start of the second data assimilation loop.
      if (iter.eq.max_iter) then
        write(niter,'(i5.5)') 0
        write(iloop,'(i1.1)') 2
        fname = name1//niter//name2//iloop//name3

c Save the data.
        open(151,file=fname,
     &    form='unformatted',access='direct',recl=4*nx*ny)

        write (151,rec=1) ((snow_d(i,j),i=1,nx),j=1,ny)
        write (151,rec=2) ((snow_depth(i,j),i=1,nx),j=1,ny)
        write (151,rec=3) ((canopy_int(i,j),i=1,nx),j=1,ny)
        write (151,rec=4) ((soft_snow_d(i,j),i=1,nx),j=1,ny)
        write (151,rec=5) ((ro_snow_grid(i,j),i=1,nx),j=1,ny)
        write (151,rec=6) ((swe_depth(i,j),i=1,nx),j=1,ny)
        write (151,rec=7) ((ro_soft_snow_old(i,j),i=1,nx),j=1,ny)
        write (151,rec=8) ((snow_d_init(i,j),i=1,nx),j=1,ny)
        write (151,rec=9) ((swe_depth_old(i,j),i=1,nx),j=1,ny)
        write (151,rec=10) ((canopy_int_old(i,j),i=1,nx),j=1,ny)
        write (151,rec=11) ((topo(i,j),i=1,nx),j=1,ny)
        write (151,rec=12) ((sum_sprec(i,j),i=1,nx),j=1,ny)

        close (151)
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine HRESTART_SAVE_DA(nx,ny,max_iter,corr_factor,
     &  icorr_factor_index,nobs_dates)

      implicit none

      include 'snowmodel.inc'

      integer iobs_num,nobs_dates,nx,ny,i,j,max_iter,iter
      integer icorr_factor_index(max_time_steps)
      real corr_factor(nx_max,ny_max,max_obs_dates+1)


c Save the correction factors for each observation date.
      open(152,file='hrestart/hrestart_corr_factor.gdat',
     &  form='unformatted',access='direct',recl=4*nx*ny)

      do iobs_num=1,nobs_dates+1
        write(152,rec=iobs_num)
     &    ((corr_factor(i,j,iobs_num),i=1,nx),j=1,ny)
      enddo

      close (152)

c Save the correction factor index.
      open(153,file='hrestart/hrestart_corr_factor_index.dat')

      do iter=1,max_iter
        write (153,*) iter,icorr_factor_index(iter)
      enddo

      close (153)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine HRESTART_READ(nx,ny,snow_d,snow_depth,
     &  canopy_int,soft_snow_d,ro_snow_grid,swe_depth,
     &  ro_soft_snow_old,snow_d_init,swe_depth_old,
     &  canopy_int_old,topo,sum_sprec,ihrestart_flag,
     &  i_dataassim_loop)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,ihrestart_flag,i_dataassim_loop,
     &  i_dataassim_loop_tmp
      real snow_d(nx_max,ny_max)
      real snow_depth(nx_max,ny_max)
      real canopy_int(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)
      real ro_snow_grid(nx_max,ny_max)
      real swe_depth(nx_max,ny_max)
      real ro_soft_snow_old(nx_max,ny_max)
      real snow_d_init(nx_max,ny_max)
      real swe_depth_old(nx_max,ny_max)
      real canopy_int_old(nx_max,ny_max)
      real topo(nx_max,ny_max)
      real sum_sprec(nx_max,ny_max)

      character*18 name1
      character*1 name2
      character*5 name3
      character*5 niter
      character*1 iloop
      character*30 fname

c Build the file name so it includes the interation number.
      name1 = 'hrestart/hrestart_'
      name2 = '_'
      name3 = '.gdat'

      if (i_dataassim_loop.lt.0.0) then
        i_dataassim_loop_tmp = 2
      else
        i_dataassim_loop_tmp = 1
      endif

      write(niter,'(i5.5)') ihrestart_flag
      write(iloop,'(i1.1)') i_dataassim_loop_tmp
      fname = name1//niter//name2//iloop//name3

c Save the data.
      open(152,file=fname,
     &  form='unformatted',access='direct',recl=4*nx*ny)

      read (152,rec=1) ((snow_d(i,j),i=1,nx),j=1,ny)
      read (152,rec=2) ((snow_depth(i,j),i=1,nx),j=1,ny)
      read (152,rec=3) ((canopy_int(i,j),i=1,nx),j=1,ny)
      read (152,rec=4) ((soft_snow_d(i,j),i=1,nx),j=1,ny)
      read (152,rec=5) ((ro_snow_grid(i,j),i=1,nx),j=1,ny)
      read (152,rec=6) ((swe_depth(i,j),i=1,nx),j=1,ny)
      read (152,rec=7) ((ro_soft_snow_old(i,j),i=1,nx),j=1,ny)
      read (152,rec=8) ((snow_d_init(i,j),i=1,nx),j=1,ny)
      read (152,rec=9) ((swe_depth_old(i,j),i=1,nx),j=1,ny)
      read (152,rec=10) ((canopy_int_old(i,j),i=1,nx),j=1,ny)
      read (152,rec=11) ((topo(i,j),i=1,nx),j=1,ny)
      read (152,rec=12) ((sum_sprec(i,j),i=1,nx),j=1,ny)

      close (152)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine HRESTART_READ_DA(nx,ny,max_iter,corr_factor,
     &  icorr_factor_index,i_dataassim_loop)

      implicit none

      include 'snowmodel.inc'

      integer iobs_num,nobs_dates,nx,ny,i,j,max_iter,iter,dummy,
     &  i_dataassim_loop
      integer icorr_factor_index(max_time_steps)
      real corr_factor(nx_max,ny_max,max_obs_dates+1)

c Read the correction factors for each observation date.
      open(152,file='hrestart/hrestart_corr_factor.gdat',
     &  form='unformatted',access='direct',recl=4*nx*ny)

      if (i_dataassim_loop.lt.0.0) nobs_dates = - i_dataassim_loop

      do iobs_num=1,nobs_dates+1
        read(152,rec=iobs_num)
     &    ((corr_factor(i,j,iobs_num),i=1,nx),j=1,ny)
      enddo

      close (152)

c Read the correction factor index.
      open(153,file='hrestart/hrestart_corr_factor_index.dat')

      do iter=1,max_iter
        read (153,*) dummy,icorr_factor_index(iter)
      enddo

      close (153)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

