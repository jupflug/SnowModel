ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c snowmodel_vars.inc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c SnowTran-3D variables.
      integer max_iter

      integer i,j,iter,nx,ny

      real ro_snow,ro_water,ro_air,gravity,vonKarman,snow_z0
      real deltax,deltay,dt,curvewt,utau_t_flag
      real fetch,xmu,C_z,h_const,wind_min,windspd_flag
      real Up_const,dz_susp,ztop_susp,fall_vel,Ur_const
      real pi,twopio360,bc_flag,topoflag,Utau_t_const
      real ht_windobs,ht_rhobs,slopewt,bs_flag,twolayer_flag
      real subgrid_flag,erosion_dist,tp_scale,curve_len_scale
      real snow_d_init_const,const_veg_flag,winddir_flag
      real windspd_min,tabler_dir,slope_adjust

      real topo_land(nx_max,ny_max)
      real topo_land2(nx_max,ny_max)
      real tabler_nn(nx_max,ny_max)
      real tabler_ss(nx_max,ny_max)
      real tabler_ee(nx_max,ny_max)
      real tabler_ww(nx_max,ny_max)
      real tabler_ne(nx_max,ny_max)
      real tabler_se(nx_max,ny_max)
      real tabler_sw(nx_max,ny_max)
      real tabler_nw(nx_max,ny_max)
      real topo(nx_max,ny_max)
      real vegtype(nx_max,ny_max)

      real uwind_grid(nx_max,ny_max),vwind_grid(nx_max,ny_max)
      real windspd_grid(nx_max,ny_max),winddir_grid(nx_max,ny_max)
      real tair_grid(nx_max,ny_max),sprec(nx_max,ny_max)
      real rh_grid(nx_max,ny_max),curvature(nx_max,ny_max)

      real slope_az(nx_max,ny_max)
      real terrain_slope(nx_max,ny_max)

      integer index_ue(ny_max,2*nx_max+1),index_uw(ny_max,2*nx_max+1)
      integer index_vn(nx_max,2*ny_max+1),index_vs(nx_max,2*ny_max+1)

      real snow_d(nx_max,ny_max)
      real snow_d_init(nx_max,ny_max)
      real Utau(nx_max,ny_max)
      real Utau_t(nx_max,ny_max)
      real z_0(nx_max,ny_max)
      real h_star(nx_max,ny_max)
      real conc_salt(nx_max,ny_max)

      real Qsalt_max(nx_max,ny_max)
      real Qsalt_maxu(nx_max,ny_max),Qsalt_maxv(nx_max,ny_max)
      real Qsalt(nx_max,ny_max)
      real Qsalt_u(nx_max,ny_max),Qsalt_v(nx_max,ny_max)
      real dh_salt(nx_max,ny_max)
      real dh_salt_u(nx_max,ny_max),dh_salt_v(nx_max,ny_max)

      real Qsusp(nx_max,ny_max)
      real Qsusp_u(nx_max,ny_max),Qsusp_v(nx_max,ny_max)
      real dh_susp(nx_max,ny_max)
      real dh_susp_u(nx_max,ny_max),dh_susp_v(nx_max,ny_max)

      real Qsubl(nx_max,ny_max)

      real sum_sprec(nx_max,ny_max)
      real wbal_qsubl(nx_max,ny_max)
      real wbal_salt(nx_max,ny_max)
      real wbal_susp(nx_max,ny_max)
      real wbal_subgrid(nx_max,ny_max)
      real sum_qsubl(nx_max,ny_max)
      real sum_trans(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)
      real ro_soft_snow(nx_max,ny_max)
      real ro_soft_snow_old(nx_max,ny_max)

      real vegsnowdepth(nvegtypes)
      real veg_z0(nx_max,ny_max)
      real vegsnowd_xy(nx_max,ny_max)
      integer iveg_ht_flag

      real run_micromet,run_enbal,run_snowpack,run_snowtran
      real print_micromet,print_enbal,print_snowpack,print_snowtran
      real print_user,print_inc
      integer iprint_inc

      character*80 topoveg_fname,met_input_fname,topo_ascii_fname,
     &  veg_ascii_fname

      character*80 snowtran_output_fname

      real curve_lg_scale_flag
      real curve_wt_lg(nx_max,ny_max)
c End SnowTran-3D variables.

c MicroMet variables.
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell
      double precision xg_line(nx_max,ny_max),yg_line(nx_max,ny_max)

      real dn                  ! average observation spacing

      real Qsi_grid(nx_max,ny_max)    ! output
      real Qli_grid(nx_max,ny_max)    ! output
      real prec_grid(nx_max,ny_max)   ! output
      real prec_grid_sol(nx_max,ny_max)
      real xlat_grid(nx_max,ny_max)   ! lat (dec deg) of cell centers
      real xlon_grid(nx_max,ny_max)   ! lon (dec deg) of cell centers

      integer iyear_init     ! model start year
      integer imonth_init    ! model start month
      integer iday_init      ! model start day
      real xhour_init        ! model start hour
      real xlat      ! approx. latitude of domain center, decimal deg

      real undef       ! undefined value
      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer i_tair_flag,i_rh_flag,i_wind_flag,i_solar_flag,
     &  i_prec_flag,i_longwave_flag,isingle_stn_flag,igrads_metfile,
     &  lapse_rate_user_flag,iprecip_lapse_rate_user_flag,n_stns_used,
     &  lat_solar_flag,ihrestart_inc,iter_start,ihrestart_flag,
     &  albedo_flag

      character*80 micromet_output_fname
      real xhour,ascii_topoveg,use_shortwave_obs,gap_frac,
     &  use_longwave_obs,use_sfc_pressure_obs,calc_subcanopy_met,
     &  cloud_frac_factor,barnes_lg_domain,UTC_flag,check_met_data,
     &  snowmodel_line_flag,wind_lapse_rate

      integer nftypes
      parameter (nftypes=5)
      real forest_LAI(nftypes)

      integer iyear,imonth,iday
      integer k_stn(nx_max,ny_max,5)
c End MicroMet variables.

c EnBal variables.
      integer icond_flag

      real Tsfc(nx_max,ny_max),Qle(nx_max,ny_max),
     &  Qh(nx_max,ny_max),Qe(nx_max,ny_max),Qc(nx_max,ny_max),
     &  Qm(nx_max,ny_max),e_balance(nx_max,ny_max),Qf(nx_max,ny_max),
     &  swe_depth(nx_max,ny_max),sfc_pressure(nx_max,ny_max),
     &  albedo(nx_max,ny_max),Qp(nx_max,ny_max)

      real albedo_glacier,albedo_diff,al_max,al_min,al_dec_cold,
     &  al_dec_melt,emiss_sfc

      character*80 enbal_output_fname
c End EnBal variables.

c SnowPack variables.
      real ro_nsnow(nx_max,ny_max),w_balance(nx_max,ny_max),
     &  runoff(nx_max,ny_max),rain(nx_max,ny_max),
     &  sum_prec(nx_max,ny_max),sum_runoff(nx_max,ny_max),
     &  xro_snow(nx_max,ny_max),ro_snow_grid(nx_max,ny_max),
     &  sum_Qcs(nx_max,ny_max),canopy_int(nx_max,ny_max),
     &  Qcs(nx_max,ny_max),canopy_unload(nx_max,ny_max),
     &  snow_depth(nx_max,ny_max),glacier_melt(nx_max,ny_max),
     &  sum_unload(nx_max,ny_max),sum_glacmelt(nx_max,ny_max),
     &  swemelt(nx_max,ny_max),d_canopy_int(nx_max,ny_max),
     &  sum_d_canopy_int(nx_max,ny_max),sum_sfcsublim(nx_max,ny_max),
     &  sum_swemelt(nx_max,ny_max),swesublim(nx_max,ny_max),
     &  swe_depth_old(nx_max,ny_max),canopy_int_old(nx_max,ny_max),
     &  liqfrac(nx_max,ny_max),icefrac(nx_max,ny_max),
     &  Tsnow(nx_max,ny_max),multResid(nx_max,ny_max),
     &  pore_space(nx_max,ny_max),mLayerThetaResid(nx_max,ny_max),
     &  availCap(nx_max,ny_max),relSaturn(nx_max,ny_max),
     &  mLayerVolFracLiqTrial(nx_max,ny_max)

      real sfc_sublim_flag,fc_param,t_avg

      character*80 snowpack_output_fname,multilayer_output_fname

      integer max_layers,multilayer_snowpack,k
      integer JJ(nx_max,ny_max)
      integer melt_flag(nx_max,ny_max,nz_max)

      real ro_snowmax,tsls_threshold,dz_snow_min,print_multilayer
      real tslsnowfall(nx_max,ny_max)
      real change_layer(nx_max,ny_max)
      real dy_snow(nx_max,ny_max,nz_max)
      real swe_lyr(nx_max,ny_max,nz_max)
      real ro_layer(nx_max,ny_max,nz_max)
      real T_old(nx_max,ny_max,nz_max)
      real gamma(nx_max,ny_max,nz_max)
      real layerFlux(nx_max,ny_max,nz_max)

      integer iclear_mn,iclear_dy,izero_snow_date
      real xclear_hr
c End SnowPack variables.

c SeaIce variables.
      real seaice_run
      real seaice_conc(nx_max,ny_max,nz_max)
c End SeaIce variables.

c Data assimilaion (precipitation and melt) factor variables.
      real corr_factor(nx_max,ny_max,max_obs_dates+1)
      integer icorr_factor_index(max_time_steps)

      integer icorr_factor_loop,irun_data_assim,nobs_dates,
     & i_dataassim_loop,i_corr_start,depth_assim
c End data assimilaion (precipitation and melt) factor variables.

      real pertPrec
c SNOWSLIDE routine
      real accumlayersum(nx_max,ny_max)
      real accumlayer(nx_max,ny_max)
      real alldat(nx_max*ny_max,7)
      real storage(nx_max,ny_max)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

