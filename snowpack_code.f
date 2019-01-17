c snowpack_code.f

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SNOWPACK_CODE(nx,ny,Tair_grid,rh_grid,ro_nsnow,
     &  dt,swe_depth,Tsfc,snow_d,prec_grid,runoff,Qm,rain,
     &  sprec,iter,w_balance,sum_prec,sum_runoff,xro_snow,
     &  undef,ro_snow,ro_snow_grid,soft_snow_d,sum_sprec,
     &  snow_depth,windspd_grid,Qsi_grid,sum_Qcs,canopy_int,
     &  Qcs,vegtype,forest_LAI,albedo,glacier_melt,
     &  canopy_unload,sum_unload,sum_glacmelt,run_snowtran,
     &  swemelt,d_canopy_int,sum_d_canopy_int,snow_d_init,
     &  sfc_pressure,Qe,sfc_sublim_flag,sum_sfcsublim,
     &  sum_swemelt,corr_factor,icorr_factor_index,swesublim,
     &  swe_depth_old,canopy_int_old,JJ,max_layers,melt_flag,
     &  ro_snowmax,tsls_threshold,dz_snow_min,tslsnowfall,
     &  change_layer,dy_snow,swe_lyr,ro_layer,T_old,gamma,
     &  multilayer_snowpack,seaice_run,seaice_conc,
     &  fc_param,t_avg,Saturn)

      implicit none

      include 'snowmodel.inc'

      integer nx,ny,iter,i,j,max_iter

      integer max_layers,multilayer_snowpack,k,n_tsteps_in_day,irec
      integer JJ(nx_max,ny_max)
      integer melt_flag(nx_max,ny_max,nz_max)

      real ro_snowmax,tsls_threshold,dz_snow_min,Cp_snow,
     &  fc_param,t_avg
      real tslsnowfall(nx_max,ny_max)
      real change_layer(nx_max,ny_max)
      real dy_snow(nx_max,ny_max,nz_max)
      real swe_lyr(nx_max,ny_max,nz_max)
      real ro_layer(nx_max,ny_max,nz_max)
      real T_old(nx_max,ny_max,nz_max)
      real gamma(nx_max,ny_max,nz_max)

      real layerFlux(nx_max,ny_max,nz_max)
      real ml_ret(nx_max,ny_max,nz_max)
      real liqfracml(nx_max,ny_max,nz_max)
      real icefracml(nx_max,ny_max,nz_max)
      real Saturn(nx_max,ny_max,nz_max)

      integer melt_flag_z(nz_max)
      real dy_snow_z(nz_max)
      real swe_lyr_z(nz_max)
      real ro_layer_z(nz_max)
      real T_old_z(nz_max)
      real gamma_z(nz_max)
      real ml_ret_z(nz_max)
      real Saturn_z(nz_max)

      real Tair_grid(nx_max,ny_max)
      real rh_grid(nx_max,ny_max)
      real prec_grid(nx_max,ny_max)
      real windspd_grid(nx_max,ny_max)
      real Qsi_grid(nx_max,ny_max)
      real vegtype(nx_max,ny_max)
      real albedo(nx_max,ny_max)
      real glacier_melt(nx_max,ny_max)
      real canopy_unload(nx_max,ny_max)
      real sum_unload(nx_max,ny_max)
      real sum_glacmelt(nx_max,ny_max)
      real sum_swemelt(nx_max,ny_max)
      real swemelt(nx_max,ny_max)
      real swesublim(nx_max,ny_max)
      real snow_d_init(nx_max,ny_max)
      real swe_depth_old(nx_max,ny_max)
      real canopy_int_old(nx_max,ny_max)
      real seaice_conc(nx_max,ny_max)

      real ro_nsnow(nx_max,ny_max),snow_d(nx_max,ny_max),
     &  runoff(nx_max,ny_max),rain(nx_max,ny_max),
     &  sprec(nx_max,ny_max),w_balance(nx_max,ny_max),
     &  sum_prec(nx_max,ny_max),sum_runoff(nx_max,ny_max),
     &  xro_snow(nx_max,ny_max),sfc_pressure(nx_max,ny_max),
     &  ro_snow_grid(nx_max,ny_max),swe_depth(nx_max,ny_max),
     &  Tsfc(nx_max,ny_max),Qm(nx_max,ny_max),
     &  soft_snow_d(nx_max,ny_max),sum_sprec(nx_max,ny_max),
     &  ro_snow,snow_depth(nx_max,ny_max),sum_Qcs(nx_max,ny_max),
     &  canopy_int(nx_max,ny_max),Qcs(nx_max,ny_max),
     &  d_canopy_int(nx_max,ny_max),sum_d_canopy_int(nx_max,ny_max),
     &  Qe(nx_max,ny_max),sum_sfcsublim(nx_max,ny_max)

      real dt,undef,Cp,xLf,Tf,A1,A2,ro_water,xLs,ro_ice,Twb,
     &  run_snowtran,sfc_sublim_flag,seaice_run,Cp_water

      real corr_factor(nx_max,ny_max,max_obs_dates+1)
      integer icorr_factor_index(max_time_steps)

      integer nftypes
      parameter (nftypes=5)
      real forest_LAI(nftypes)

      print *,'   solving the snow-cover evolution'

c Define the constants used in the computations.
      CALL CONSTS_SNOWPACK(Cp,xLs,ro_ice,xLf,Tf,A1,A2,ro_water,
     &  Cp_snow,ro_snowmax,Cp_water)

c Run the snowpack evolution sub-model.
      do j=1,ny
        do i=1,nx

c Extract the vertical column for this i,j point, and send it
c   to the subroutine. *** Note that I should use f95, then I would
c   not have to do this (I could pass in subsections of the arrays).
          if (multilayer_snowpack.eq.1) then
            do k=1,nz_max
              melt_flag_z(k) = melt_flag(i,j,k)
              dy_snow_z(k) = dy_snow(i,j,k)
              swe_lyr_z(k) = swe_lyr(i,j,k)
              ro_layer_z(k) = ro_layer(i,j,k)
              T_old_z(k) = T_old(i,j,k)
              gamma_z(k) = gamma(i,j,k)

c J.PFLUG
c variables for multilayer percolation investigation
c              layerFlux_z(k) = layerFlux(i,j,k)
              ml_ret_z(k) = ml_ret(i,j,k)
              Saturn_z(k) = Saturn(i,j,k)
c END J.PFLUG
            enddo
          endif

          CALL SNOWPACK_CORE(Twb,Tf,Tair_grid(i,j),rh_grid(i,j),xLs,
     &      Cp,sfc_pressure(i,j),ro_nsnow(i,j),dt,ro_snow,
     &      swe_depth(i,j),Tsfc(i,j),A1,A2,snow_d(i,j),ro_water,
     &      ro_ice,prec_grid(i,j),runoff(i,j),Qm(i,j),xLf,rain(i,j),
     &      sprec(i,j),iter,w_balance(i,j),sum_prec(i,j),
     &      sum_runoff(i,j),xro_snow(i,j),undef,
     &      soft_snow_d(i,j),sum_sprec(i,j),ro_snow_grid(i,j),
     &      snow_depth(i,j),windspd_grid(i,j),Qsi_grid(i,j),
     &      sum_Qcs(i,j),canopy_int(i,j),Qcs(i,j),vegtype(i,j),
     &      forest_LAI,albedo(i,j),canopy_unload(i,j),
     &      sum_unload(i,j),sum_glacmelt(i,j),run_snowtran,
     &      swemelt(i,j),d_canopy_int(i,j),sum_d_canopy_int(i,j),
     &      snow_d_init(i,j),Qe(i,j),glacier_melt(i,j),
     &      sfc_sublim_flag,sum_sfcsublim(i,j),sum_swemelt(i,j),
     &      corr_factor(i,j,-icorr_factor_index(iter)),
     &      icorr_factor_index(iter),swesublim(i,j),
     &      swe_depth_old(i,j),canopy_int_old(i,j),JJ(i,j),
     &      max_layers,melt_flag_z,ro_snowmax,tsls_threshold,
     &      dz_snow_min,tslsnowfall(i,j),change_layer(i,j),dy_snow_z,
     &      swe_lyr_z,ro_layer_z,T_old_z,gamma_z,multilayer_snowpack,
     &      Cp_snow,seaice_run,fc_param,t_avg,Cp_water,ml_ret_z,
     &      Saturn_z)

c Re-build the 3-D arrays.  See note above about using f95 to avoid this.
          if (multilayer_snowpack.eq.1) then
            do k=1,nz_max
              melt_flag(i,j,k) = melt_flag_z(k)
              dy_snow(i,j,k) = dy_snow_z(k)
              swe_lyr(i,j,k) = swe_lyr_z(k)
              ro_layer(i,j,k) = ro_layer_z(k)
              T_old(i,j,k) = T_old_z(k)
              gamma(i,j,k) = gamma_z(k)
              ml_ret(i,j,k) = ml_ret_z(k)
              Saturn(i,j,k) = Saturn_z(k)
            enddo
          endif

        enddo
      enddo

      if (run_snowtran.eq.0.0) then
        do j=1,ny
          do i=1,nx
          swe_depth_old(i,j) = swe_depth(i,j)
          canopy_int_old(i,j) = canopy_int(i,j)
          enddo
        enddo
      endif

c Read in the sea ice concentration.  These are daily data, so
c   first calculate which record in the data file this time step
c   corresponds to.
      if (seaice_run.ne.0.0) then
        n_tsteps_in_day = nint(86400.0 / dt)
        if (mod(iter-1,n_tsteps_in_day).eq.0) then
          irec = int((real(iter) - 0.5) * dt / 86400.0) + 1
          print *,'sea ice irec =',irec
          read (445,rec=irec) ((seaice_conc(i,j),i=1,nx),j=1,ny)
        endif
      endif

c If this simulation is not running SnowTran-3D, then zero out
c   the ocean grid cells that have no sea ice here.  If it is
c   running with SnowTran-3D, then do this in the SnowTran-3D
c   subroutine.
      if (run_snowtran.eq.0.0) then
        if (seaice_run.ne.0.0) then
          CALL ZERO_SEAICE_SNOW(nx,ny,snow_depth,ro_snow_grid,
     &      ro_snow,swe_depth,swe_depth_old,canopy_int_old,JJ,
     &      tslsnowfall,dy_snow,swe_lyr,ro_layer,T_old,
     &      multilayer_snowpack,tsls_threshold,seaice_conc)
        endif
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SNOWPACK_CORE(Twb,Tf,Tair,rh,xLs,
     &  Cp,sfc_pressure,ro_nsnow,dt,ro_snow,
     &  swe_depth,Tsfc,A1,A2,snow_d,ro_water,
     &  ro_ice,prec,runoff,Qm,xLf,rain,
     &  sprec,iter,w_balance,sum_prec,
     &  sum_runoff,xro_snow,undef,
     &  soft_snow_d,sum_sprec,ro_snow_grid,
     &  snow_depth,windspd,Qsi,
     &  sum_Qcs,canopy_int,Qcs,vegtype,
     &  forest_LAI,albedo,canopy_unload,
     &  sum_unload,sum_glacmelt,run_snowtran,
     &  swemelt,d_canopy_int,sum_d_canopy_int,
     &  snow_d_init,Qe,glacier_melt,
     &  sfc_sublim_flag,sum_sfcsublim,sum_swemelt,
     &  corr_factor,
     &  icorr_factor_index,swesublim,
     &  swe_depth_old,canopy_int_old,JJ,
     &  max_layers,melt_flag,ro_snowmax,tsls_threshold,
     &  dz_snow_min,tslsnowfall,change_layer,dy_snow,
     &  swe_lyr,ro_layer,T_old,gamma,multilayer_snowpack,
     &  Cp_snow,seaice_run,fc_param,t_avg,
     &  Cp_water,ml_ret,Saturn)

      implicit none

      include 'snowmodel.inc'

      integer iter,icorr_factor_index

      integer JJ,max_layers,multilayer_snowpack
      
      real ro_snowmax,tsls_threshold,dz_snow_min,tslsnowfall,Cp_snow,
     &  Cp_water

      integer melt_flag(nz_max)
      real change_layer
      real dy_snow(nz_max)
      real swe_lyr(nz_max)
      real ro_layer(nz_max)
      real T_old(nz_max)
      real gamma(nz_max)
      real ml_ret(nz_max)
      real Saturn(nz_max)

      real Twb,Tf,Tair,rh,xLs,Cp,ro_nsnow,dt,ro_snow,swe_depth,
     &  Tsfc,A1,A2,snow_d,ro_water,ro_ice,prec,runoff,Qm,xLf,rain,
     &  sprec,w_balance,sum_prec,sum_runoff,xro_snow,undef,
     &  soft_snow_d,sum_sprec,ro_snow_grid,snow_depth,sprec_grnd,
     &  windspd,Qsi,sum_Qcs,canopy_int,Qcs,canopy_unload,
     &  vegtype,albedo,glacier_melt,sum_unload,sum_glacmelt,
     &  run_snowtran,swemelt,d_canopy_int,sfc_pressure,
     &  sum_d_canopy_int,snow_d_init,Qe,sfc_sublim_flag,
     &  sum_sfcsublim,sum_swemelt,corr_factor,swesublim,
     &  swe_depth_old,canopy_int_old,sprec_grnd_ml,seaice_run,
     &  mLayerVolFracLiqTrial,fc_param,t_avg

      integer nftypes
      parameter (nftypes=5)
      real forest_LAI(nftypes)

c Calculate the canopy sublimation, loading and unloading.  Note
c   that here I have assumed that evergreen trees are type 1.
      if (vegtype.le.5.0) then
        CALL CANOPY_SNOW(rh,Tair,windspd,Qsi,sum_Qcs,albedo,
     &    canopy_int,sprec,Qcs,dt,canopy_unload,
     &    forest_LAI(nint(vegtype)),sum_unload,d_canopy_int,
     &    sum_d_canopy_int)
        sprec_grnd = sprec + canopy_unload - d_canopy_int
        sprec_grnd_ml = sprec - d_canopy_int
      else
        Qcs = 0.0
        sprec_grnd = sprec
        sprec_grnd_ml = sprec
      endif

c Solve for the wet bulb temperature.
      CALL SOLVEWB(Twb,Tf,Tair,rh,xLs,Cp,sfc_pressure)

c Compute the new snow density.
      CALL NSNOWDEN(ro_nsnow,Twb,Tf,dt)
     
c Call the multi-layer snowpack model.
      if (multilayer_snowpack.eq.1) then

        CALL MULTI_LAYER_SNOW(JJ,ro_layer,Tf,dt,ro_water,
     &    ro_ice,T_old,dy_snow,swe_lyr,Qm,ro_snowmax,rain,
     &    xLf,Cp_snow,melt_flag,runoff,tslsnowfall,ro_nsnow,
     &    sprec,Tsfc,tsls_threshold,gamma,max_layers,change_layer,
     &    dz_snow_min,snow_depth,swe_depth,undef,canopy_unload,
     &    vegtype,glacier_melt,sum_glacmelt,sum_swemelt,snow_d,
     &    Qe,sfc_sublim_flag,sum_sfcsublim,soft_snow_d,ro_snow,
     &    sum_sprec,sprec_grnd_ml,sum_prec,prec,sum_runoff,
     &    ro_snow_grid,xro_snow,swesublim,A1,A2,
     &    fc_param,Cp_water,ml_ret,Saturn)

c Call the original single-layer snowpack model.
      else

c Compute the snow density change due to settling.
        CALL DDENSITY(ro_snow_grid,swe_depth,Tf,Tsfc,dt,A1,A2,
     &    snow_depth,ro_water,ro_ice)

c Compute the melt, rain, and snow contributions to modifying
c   the snowpack depth, density, and snow water equivalent.
        CALL SNOWPACK(swe_depth,snow_d,ro_snow_grid,
     &    prec,ro_water,ro_nsnow,runoff,Qm,xLf,dt,rain,sprec,
     &    sum_prec,sum_runoff,soft_snow_d,sum_sprec,ro_snow,
     &    snow_depth,sprec_grnd,vegtype,glacier_melt,sum_glacmelt,
     &    swemelt,canopy_unload,Qe,sfc_sublim_flag,sum_sfcsublim,
     &    sum_swemelt,corr_factor,icorr_factor_index,swesublim,
     &    ro_snowmax)

c Post process the data for output.
        CALL POSTPROC(ro_snow_grid,xro_snow,snow_depth,undef)
      endif

c Perform a water balance check (see notes in this subroutine).
      if (seaice_run.eq.0.0) then
        if (run_snowtran.eq.0.0) then
          CALL WATERBAL_SNOWPACK(w_balance,prec,Qcs,runoff,
     &    d_canopy_int,swe_depth,glacier_melt,swe_depth_old,iter,
     &    swesublim,canopy_unload,canopy_int_old,canopy_int)

c         CALL WATERBAL_SNOWPACK_sums(w_balance,sum_prec,sum_Qcs,
c    &      sum_runoff,canopy_int,swe_depth,sum_glacmelt,iter,
c    &      snow_d_init,ro_snow,ro_water,sum_sfcsublim)
        endif
      endif
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CANOPY_SNOW(rh,Tair,windspd,Qsi,sum_Qcs,albedo,
     &  canopy_int,sprec,Qcs,dt,canopy_unload,
     &  forest_LAI,sum_unload,d_canopy_int,
     &  sum_d_canopy_int)

      implicit none

      real rh,Tair,windspd,V_s,Qsi,forest_LAI,dt,xImax,canopy_int,
     &  d_canopy_int,Qcs,Ce,sprec,C_0,unload_melt,canopy_unload,
     &  sum_Qcs,albedo,sum_unload,sum_d_canopy_int

c Note that all of this must deal with the (kg/m2)=(mm), => (m)
c   issues.  Precip is in (m), all of these equations are in
c   (kg/m2), and I want the outputs to be in (m).

c Compute the sublimation loss rate coefficient for canopy snow.
      CALL SUBLIM_COEF(rh,Tair,windspd,V_s,Qsi,albedo)

c Maximum interception storage.
      xImax = 4.4 * forest_LAI

c Change in canopy load due to snow precipitation during this time
c   step.  Convert the canopy interception to mm.
      canopy_int = 1000.0 * canopy_int
      d_canopy_int = 0.7 * (xImax - canopy_int) *
     &  ( 1.0 - exp((- sprec)*1000.0/xImax))

c Update the interception load.
      canopy_int = canopy_int + d_canopy_int

c Canopy exposure coefficient.
      if (canopy_int.eq.0.0) then
        Ce = 0.0
      else
c Pomeroy's k_c value
c       Ce = 0.0114 * (canopy_int/xImax)**(-0.4)
c My k_c value.
        Ce = 0.00995 * (canopy_int/xImax)**(-0.4)
      endif

c Canopy sublimation (kg/m2), (a negative mumber).  Make sure that
c   you don't sublimate more than is available.
      Qcs = Ce * canopy_int * V_s * dt
      Qcs = -min(canopy_int,-Qcs)

c Remove the sublimated moisture from the canopy store.
      canopy_int = canopy_int + Qcs

c Save the sublimation in (m).
      Qcs = Qcs / 1000.0
      sum_Qcs = sum_Qcs + Qcs

c Perform a second unloading due to melt.  Assume an unloading rate
c   of 5.0 mm/day/C.
      C_0 = 5.0 / 86400.0
      unload_melt = C_0 * max(0.0,Tair-273.16) * dt
      unload_melt = min(canopy_int,unload_melt)
      canopy_int = canopy_int - unload_melt

c Keep track of the unloaded snow that reached the ground during
c   this time step (m) (this will add to the snow depth).
      canopy_unload = unload_melt / 1000.0
      d_canopy_int = d_canopy_int / 1000.0

c Save a summing array of this unloaded snow.
      sum_unload = sum_unload + canopy_unload

c Save a summing array of the change in canopy load.
      sum_d_canopy_int = sum_d_canopy_int + d_canopy_int

c Save the interception load for the next time step.  Convert to m.
      canopy_int = canopy_int / 1000.0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SUBLIM_COEF(rh,Tair,windspd,V_s,Qsi,albedo)

c Compute the sublimation loss rate coefficient for canopy snow.

      implicit none

      real pi,ro_ice,xM,R,R_dryair,vonKarman,visc_air,h_s,xlamdaT,
     &  D,ro_sat,sigma,V_s,radius,xmass,windspd,rh,Tair,Qsi,Sp,
     &  xN_r,xNu,xSh,top,bottom,omega,albedo

c Constants.
      pi = 2.0 * acos(0.0)
      ro_ice = 917.0
      xM = 18.01
      R = 8313.
      R_dryair = 287.
      vonKarman = 0.4
      visc_air = 13.e-6
      h_s = 2.838e6
      xlamdaT = 0.024

c Particle radius.
      radius = 5.0e-4

c Particle mass.
      xmass = 4.0/3.0 * pi * ro_ice * radius**3

c Diffusivity of water vapor in the atmosphere.
      D = 2.06e-5 * (Tair/273.)**(1.75)

c Saturation density of water vapor.
      ro_sat = 0.622 / (R_dryair * Tair) *
     &  611.15 * exp(22.452 * (Tair - 273.16) / (Tair - 0.61))

c Humidity deficit.
      sigma = 0.01 * rh - 1.0
      sigma = min(0.0,sigma)
      sigma = max(-1.0,sigma)

c Reynolds, Nusselt, and Sherwood numbers.
      xN_r = 2.0 * radius * windspd / visc_air
      xNu = 1.79 + 0.606 * xN_r**(0.5)
      xSh = xNu

c Solar radiation absorbed by the snow particle.  Here assume that
c   the general snow albedo is the same as the snow particle albedo.
      Sp = pi * radius**2 * (1.0 - albedo) * Qsi

c Sublimation-loss rate coefficient for an ice sphere.
      omega = ((h_s * xM)/(R * Tair) - 1.0) / (xlamdaT * Tair * xNu)
      top = 2.0 * pi * radius * sigma - Sp * omega
      bottom = h_s * omega + 1.0/(D * ro_sat * xSh)
      V_s = (top/bottom)/xmass

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WATERBAL_SNOWPACK(w_balance,prec,Qcs,runoff,
     &  d_canopy_int,swe_depth,glacier_melt,swe_depth_old,iter,
     &  swesublim,canopy_unload,canopy_int_old,canopy_int)

      implicit none

      integer iter

      real w_balance,prec,Qcs,runoff,d_canopy_int,swe_depth_old,
     &  swe_depth,glacier_melt,swesublim,canopy_unload,canopy_int_old,
     &  canopy_int

c Note that the following balances should hold.  These aren't quite
c   right, but it is a place to start.
c   Canopy Balance (forest):
c     canopy = sprec - unload + Qcs ==> unload = sprec - canopy + Qcs
c
c   Snowpack Balance (forest):
c     swe_d = unload + rain - runoff ==>
c       canopy + swe_d = sprec + rain + Qcs - runoff
c     prec = sprec + rain
c     sum_rain  = sum_sprec - sum_prec
c
c   Snowpack Balance (non-forest):
c     swe_d = sprec + rain - runoff + subl + salt + susp + subgrid +
c       glaciermelt
c
c   Everywhere:
c     w_balance = sum_prec + sum_Qcs - sum_runoff + sum_subl +
c       sum_trans - canopy_int - swe_depth + sum_glacmelt
c
c   The related variables that would need to be brought in are:
c      d_canopy_int,sum_d_canopy_int,sum_unload

c This subroutine is called for the case where SnowTran-3D is not
c   run.  The subroutine WATERBAL_SNOWTRAN is used if the model
c   simulation includes SnowTran-3D.
c     w_balance = swe_depth_old - swe_depth + prec - runoff +
c    &  glacier_melt - swesublim + canopy_int_old - canopy_int -
c    &  d_canopy_int + Qcs + canopy_unload

c Do the snowpack.
c     w_balance = swe_depth_old - swe_depth + prec - runoff -
c    &  glacier_melt - swesublim

c Do the canopy.
c     w_balance = canopy_int_old - canopy_int + d_canopy_int +
c    &  Qcs - canopy_unload

c Do the snowpack and canopy store.
      w_balance = swe_depth_old - swe_depth + prec - runoff +
     &  glacier_melt - swesublim + canopy_int_old - canopy_int +
     &  Qcs

      if (abs(w_balance).gt.1.0e-5) then
        print*,'water imbalance found, iter =',iter,' ',w_balance
        print*,swe_depth_old,swe_depth,prec,runoff,glacier_melt, 
     &    swesublim,canopy_int_old,canopy_int,Qcs
c        stop
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WATERBAL_SNOWPACK_sums(w_balance,sum_prec,sum_Qcs,
     &  sum_runoff,canopy_int,swe_depth,sum_glacmelt,iter,
     &  snow_d_init,ro_snow,ro_water,sum_sfcsublim)

      implicit none

      integer iter

      real w_balance,sum_prec,sum_Qcs,sum_runoff,canopy_int,
     &  swe_depth,sum_glacmelt,snow_d_init,ro_snow,ro_water,
     &  sum_sfcsublim

c Note that the following balances should hold.  These aren't quite
c   right, but it is a place to start.
c   Canopy Balance (forest):
c     canopy = sprec - unload + Qcs ==> unload = sprec - canopy + Qcs
c
c   Snowpack Balance (forest):
c     swe_d = unload + rain - runoff ==>
c       canopy + swe_d = sprec + rain + Qcs - runoff
c     prec = sprec + rain
c     sum_rain  = sum_sprec - sum_prec
c
c   Snowpack Balance (non-forest):
c     swe_d = sprec + rain - runoff + subl + salt + susp + subgrid +
c       glaciermelt
c
c   Everywhere:
c     w_balance = sum_prec + sum_Qcs - sum_runoff + sum_subl +
c       sum_trans - canopy_int - swe_depth + sum_glacmelt
c
c   The related variables that would need to be brought in are:
c      d_canopy_int,sum_d_canopy_int,sum_unload

c This subroutine is called for the case where SnowTran-3D is not
c   run.  The subroutine WATERBAL_SNOWTRAN is used if the model
c   simulation includes SnowTran-3D.
      w_balance = sum_prec + sum_Qcs - sum_runoff - canopy_int -
     &  swe_depth + sum_glacmelt + snow_d_init * ro_snow/ro_water -
     &  sum_sfcsublim

      if (abs(w_balance).gt.1.0e-4)
     &  print*,'water imbalance found, iter =',iter,' ',w_balance

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc;cccccc

      SUBROUTINE SNOWPACK(swe_depth,snow_d,ro_snow_grid,
     &  prec,ro_water,ro_nsnow,runoff,Qm,xLf,dt,rain,sprec,
     &  sum_prec,sum_runoff,soft_snow_d,sum_sprec,ro_snow,
     &  snow_depth,sprec_grnd,vegtype,glacier_melt,sum_glacmelt,
     &  swemelt,canopy_unload,Qe,sfc_sublim_flag,sum_sfcsublim,
     &  sum_swemelt,corr_factor,icorr_factor_index,swesublim,
     &  ro_snowmax)

      implicit none

      real ro_snowmax,runoff,Qm,swe_depth,potmelt,swemelt,dt,
     &  ro_water,xLf,snow_depth,ro_snow_grid,snow_d_melt,dz_water,
     &  soft_snow_d,prec,rain,snow_d,sum_sprec,sum_prec,
     &  sum_runoff,ro_nsnow,sprec,ro_snow,snow_d_new,sprec_grnd,
     &  vegtype,glacier_melt,sum_glacmelt,canopy_unload,Qe,
     &  xLsublim,potsublim,swesublim,snow_d_sublim,sfc_sublim_flag,
     &  sum_sfcsublim,sum_swemelt,corr_factor,potmelt_tmp
      integer icorr_factor_index

      runoff = 0.0

c SURFACE SUBLIMATION.

c Whether static-surface (non-blowing snow) sublimation is included
c   in the model calculations is controlled by the sfc_sublim_flag.
c   I am waiting for the flux-tower data Matthew and I are collecting
c   in Alaska, to compare with the model simulations, before
c   including this part of the model in all simulations.

c If the sfc_sublim_flag is turned on, the latent heat flux (Qe)
c   calculated in ENBAL is used to add/remove snow from the snowpack.
c   xLsublim = xLf + xLv = 2.5104x10^6 J/kg + 3.334x10^5 J/kg, and
c   potsublim is in m swe.

      if (swe_depth.gt.0.0  .and.  sfc_sublim_flag.eq.1.0) then
        if (Qe.lt.0.0) then

c Compute the snow-surface sublimation (m, swe).
          xLsublim = 2.844e6
          potsublim = (- dt) * Qe / (ro_water * xLsublim)
          swesublim = min(potsublim,swe_depth)

c Save a summing array of the static surface snow sublimation.
          sum_sfcsublim = sum_sfcsublim + swesublim

c Compute the change in snow depth.  Assume that this sublimated
c   snow does not change the snow density and does not change the
c   soft snow depth.  It only reduces the snow depth and the
c   associated swe depth.
          swe_depth = swe_depth - swesublim
          if (swe_depth.eq.0.0) then
            snow_depth = 0.0
          else
            snow_d_sublim = swesublim * ro_water / ro_snow_grid
            snow_depth = snow_depth - snow_d_sublim
          endif
        else
          swesublim = 0.0
        endif
      else
        swesublim = 0.0
      endif

c MELTING.

c If melting occurs, decrease the snow depth, and place the melt
c   water in the 'runoff' variable.  Keep track of the liquid water
c   produced.

      if (Qm.gt.0.0) then

c Compute the snow melt (m).
        potmelt = dt * Qm / (ro_water * xLf)

c Account for any snowmelt data assimilation.
        if (icorr_factor_index.lt.0.0) then
          potmelt_tmp = potmelt * corr_factor
          swemelt = min(potmelt_tmp,swe_depth)
c Handle the case of no snowmelt data assimilation.
        else
          swemelt = min(potmelt,swe_depth)
        endif

c Compute any glacier or permanent snow-field melt (m water equiv.).
        if (vegtype.eq.20.0) then
          glacier_melt = potmelt - swemelt
        else
          glacier_melt = 0.0
        endif

c Save a summing array of the glacier melt.
        sum_glacmelt = sum_glacmelt + glacier_melt

c Save the runoff contribution.
        runoff = runoff + glacier_melt

c Save a summing array of the snow melt.
        sum_swemelt = sum_swemelt + swemelt

c Compute the change in snow depth.
        snow_d_melt = swemelt * ro_water / ro_snow_grid
        snow_depth = snow_depth - snow_d_melt
        snow_depth = max(0.0,snow_depth)

c Compute the changes in snow density resulting from the melt.
c   Assume that the melted snow is redistributed through the new
c   snow depth up to a maximum density.  Any additional melt water
c   is added to the runoff.
        if (snow_depth.eq.0.0) then
          ro_snow_grid = ro_snowmax
          runoff = runoff + swemelt
        else
          ro_snow_grid = swe_depth * ro_water / snow_depth
        endif

        if (ro_snow_grid.gt.ro_snowmax) then
          dz_water = snow_depth *
     &      (ro_snow_grid - ro_snowmax) / ro_water
          ro_snow_grid = ro_snowmax
          swe_depth = snow_depth * ro_snow_grid / ro_water
          runoff = runoff + dz_water
        else
          swe_depth = snow_depth * ro_snow_grid / ro_water
        endif

        soft_snow_d = 0.0

      else

c These prevent values from the previous time step from being
c   carried through to the next time step.
        swemelt = 0.0
        glacier_melt = 0.0

      endif

c PRECIPITATION.

c Precipitation falling as rain on snow contributes to a snow
c   density increase, precipitation falling as snow adds to the
c   snow depth, and rain falling on bare ground contributes to the
c   runoff.

c We have precipitation.
      if (prec.gt.0.0) then
        rain = prec - sprec

c We have rain.
        if (rain.gt.0.0) then

c Rain on snow.  Note that we can also have snow unloading here.
c   Assume this unloading is wet as rain.
          if (snow_depth.gt.0.0) then
            swe_depth = swe_depth + rain + canopy_unload
            ro_snow_grid = swe_depth * ro_water / snow_depth
            if (ro_snow_grid.gt.ro_snowmax) then
              dz_water = snow_depth * (ro_snow_grid - ro_snowmax) /
     &          ro_water
              ro_snow_grid = ro_snowmax
              swe_depth = snow_depth * ro_snow_grid / ro_water
              runoff = runoff + dz_water
            endif

c Rain on bare ground.  Assume any unloading is as wet as rain.
          else
            runoff = runoff + rain + canopy_unload
          endif

c We have snow precipitation (on either snow or bare ground).
        else
          swe_depth = swe_depth + sprec_grnd
          snow_d_new = ro_water / ro_nsnow * sprec_grnd
          snow_depth = snow_depth + snow_d_new
          ro_snow_grid = ro_water * swe_depth / snow_depth
        endif

c Here we handle the case where there is no precipitation, but
c   there is snow falling from the canopy to the snowpack.
      else
        rain = 0.0
        if (sprec_grnd.gt.0.0) then
          swe_depth = swe_depth + sprec_grnd
          snow_d_new = ro_water / ro_snow * sprec_grnd
          snow_depth = snow_depth + snow_d_new
          ro_snow_grid = ro_water * swe_depth / snow_depth
        endif
      endif

c The following are set up to be compatible with SnowTran-3D, and
c   are in snow-depth units.  The sum_sprec corrections are done
c   in the SnowTran-3D code.
      soft_snow_d = soft_snow_d + sprec_grnd * ro_water / ro_snow
      snow_d = swe_depth * ro_water / ro_snow
c     sum_sprec = sum_sprec + sprec_grnd * ro_water / ro_snow
      sum_sprec = sum_sprec + sprec_grnd

c The following are in swe-depth units.
      sum_prec = sum_prec + prec
      sum_runoff = sum_runoff + runoff

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE DDENSITY(ro_snow_grid,swe_depth,Tf,Tsfc,dt,A1,A2,
     &  snow_depth,ro_water,ro_ice)

      implicit none

      real snow_depth,Tsg,Tf,Tsnow,Tsfc,ro_snow_grid,dt,A1,A2,
     &  swe_depth_star,ro_ice,ro_water,swe_depth

      if (snow_depth.gt.0.0) then

c Assume that the snow-ground interface temperature is -1.0 C.
        Tsg = Tf - 1.0
        Tsnow = 0.5 * (Tsg + Tsfc)
        swe_depth_star= 0.5 * swe_depth
        ro_snow_grid = ro_snow_grid + dt *
     &    (A1 * swe_depth_star * ro_snow_grid *
     &    exp((- 0.08)*(Tf-Tsnow)) * exp((- A2)*ro_snow_grid))
        ro_snow_grid = min(ro_ice,ro_snow_grid)
        snow_depth = ro_water * swe_depth / ro_snow_grid

      endif

      return
      end
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SOLVEWB(xnew,Tf,Tair,rh,xLs,Cp,sfc_pressure)

      implicit none

      real A,B,C,ea,rh,Tair,Tf,tol,old,fprime,xLs,Cp,funct,xnew,
     &  sfc_pressure

      integer maxiter,i

c Coeffs for saturation vapor pressure over water (Buck 1981).
c   Note: temperatures for Buck`s equations are in deg C, and
c   vapor pressures are in mb.  Do the adjustments so that the
c   calculations are done with temperatures in K, and vapor
c   pressures in Pa.

c Over water.
        A = 6.1121 * 100.0
        B = 17.502
        C = 240.97
c Over ice.
c       A = 6.1115 * 100.0
c       B = 22.452
c       C = 272.55

c Atmospheric vapor pressure from relative humidity data.
      ea = rh / 100.0 * A * exp((B * (Tair - Tf))/(C + (Tair - Tf)))

c Solve for the wet bulb temperature.
      tol = 1.0e-2
      maxiter = 20
      old = Tair

      do i=1,maxiter
        fprime = 1.0 + xLs/Cp * 0.622/sfc_pressure * log(10.0) *
     &    2353. * (10.0**(11.40 - 2353./old)) / old**2
        funct = old - Tair + xLs/Cp * 0.622/sfc_pressure *
     &    (10.0**(11.40-2353./old) - ea)
        xnew = old - funct/fprime
        if (abs(xnew - old).lt.tol) return
        old = xnew
      end do

c If the maximum iterations are exceeded, send a message and set
c   the wet bulb temperature to the air temperature.
      write (*,102)
  102 format('max iteration exceeded when solving for Twb')
      xnew = Tair

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE NSNOWDEN(ro_nsnow,Twb,Tf,dt)

      implicit none

      real Twgmax,Tf,Twb,ro_nsnow,scalefact,dt,wt

      Twgmax = Tf + 1.0
      if (Twb.ge.258.16 .and. Twb.le.Twgmax) then
        ro_nsnow = 50. + 1.7 * (Twb - 258.16)**1.5
      elseif (Twb.lt.258.16) then
        ro_nsnow = 50.0
      else
        ro_nsnow = 158.8
      endif


c For one day time steps, this equation gives a new snow density at
c   the end of the 24 hour period which is too low, by an approximate
c   factor of X.  Thus, for a daily time step, I scale the density by
c   X before returning it to the main program.

      scalefact = 1.0
      if (dt.eq.86400.0) then
        if (ro_nsnow.le.158.8) then
          wt = 1.0 + (50.0 - ro_nsnow) / 108.8
          ro_nsnow = wt * scalefact * ro_nsnow + ro_nsnow
          ro_nsnow = min(158.8,ro_nsnow)
        endif
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE POSTPROC(ro_snow_grid,xro_snow,snow_depth,undef)

      implicit none

      real snow_depth,xro_snow,undef,ro_snow_grid

      if (snow_depth.eq.0.0) then
        xro_snow = undef
      else
        xro_snow = ro_snow_grid
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CONSTS_SNOWPACK(Cp,xLs,ro_ice,xLf,Tf,A1,A2,ro_water,
     &  Cp_snow,ro_snowmax,Cp_water)

      implicit none

      real Cp,xLs,ro_ice,xLf,Tf,A1,A2,ro_water,Cp_snow,ro_snowmax,
     &  Cp_water

      Cp = 1004.
      xLs = 2.500e6
      ro_ice = 917.0
      xLf = 3.34e5
      Tf = 273.16
      A1 = 0.0013
      A2 = 0.021
      ro_water = 1000.0
      Cp_snow = 2106.
      ro_snowmax = 550.0
      Cp_water = 4180.0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MULTI_LAYER_SNOW(JJ,ro_layer,Tf,dt,ro_water,
     &  ro_ice,T_old,dy_snow,swe_lyr,Qm,ro_snowmax,rain,
     &  xLf,Cp_snow,melt_flag,runoff,tslsnowfall,ro_nsnow,
     &  sprec,Tsfc,tsls_threshold,gamma,max_layers,change_layer,
     &  dz_snow_min,snow_depth,swe_depth,undef,canopy_unload,
     &  vegtype,glacier_melt,sum_glacmelt,sum_swemelt,snow_d,
     &  Qe,sfc_sublim_flag,sum_sfcsublim,soft_snow_d,ro_snow,
     &  sum_sprec,sprec_grnd_ml,sum_prec,prec,sum_runoff,
     &  ro_snow_grid,xro_snow,swesublim,A1,A2,fc_param,
     &  Cp_water,ml_ret,Saturn)

      implicit none

      include 'snowmodel.inc'

      integer JJ,max_layers
      integer melt_flag(nz_max)
      integer i,j

      real dy_snow(nz_max)
      real swe_lyr(nz_max)
      real ro_layer(nz_max)
      real T_old(nz_max)
      real gamma(nz_max)

      real ml_ret(nz_max)
      real Saturn(nz_max)

      real Tf,dt,ro_water,ro_ice,Qm,ro_snowmax,rain,xLf,Cp_snow,
     &  runoff,tslsnowfall,ro_nsnow,sprec,Tsfc,tsls_threshold,
     &  dz_snow_min,snow_depth,swe_depth,undef,change_layer,
     &  canopy_unload,vegtype,glacier_melt,sum_glacmelt,sum_swemelt,
     &  soft_snow_d,Qe,sfc_sublim_flag,sum_sfcsublim,snow_d,
     &  ro_snow,sum_sprec,sprec_grnd_ml,sum_prec,prec,sum_runoff,
     &  ro_snow_grid,xro_snow,swesublim,A1,A2,Tk,fc_param,Cp_water

c THIS IS THE MULTI-LAYER SNOWPACK MODEL.

c Note there is some confusion with the dy - dz notation used here.
c   In the multi-layer code 'y' is the vertical coordinate.  This is
c   a hold-over from a time when my temperature solution code had
c   y going up-down.

      
c Compute the snow density change due to compaction.
      CALL DDENSITY_ML(ro_layer,Tf,dt,ro_water,ro_ice,
     &  T_old,JJ,dy_snow,A1,A2)
      
c Calculate the rainfall from prec and sprec.
      if (prec.gt.0.0) then
        rain = prec - sprec
      else
        rain = 0.0
      endif

c J.PFLUG
c included PRECIP_ML, MERGE_LAYERS_ML, SNOWTEMP_ML, and POST_PROC_ML
c within the MELT_SNOW_ML routine for convenience with the dynamic 
c timestep. This can be easily changed by looking at the code in 
c MELT_SNOW_ML 

c Distribute surface melt and rain precipitation through the snowpack.
      CALL MELT_SNOW_ML(JJ,swe_lyr,ro_water,ro_layer,Qm,dt,
     &  dy_snow,ro_snowmax,rain,xLf,Cp_snow,Tf,T_old,melt_flag,
     &  runoff,canopy_unload,swe_depth,snow_depth,vegtype,
     &  glacier_melt,sum_glacmelt,sum_swemelt,soft_snow_d,Qe,
     &  sfc_sublim_flag,sum_sfcsublim,swesublim,fc_param,
     &  Cp_water,ml_ret,
     &  change_layer,dz_snow_min,gamma,max_layers,prec,ro_nsnow,
     &  ro_snow,ro_snow_grid,snow_d,sprec_grnd_ml,sum_prec,
     &  sum_runoff,sum_sprec,tsfc,tsls_threshold,tslsnowfall,undef,
     &  xro_snow,Saturn)
     
c Account for the accumulation of snow precipitation on the snowpack.
c      CALL PRECIP_ML(JJ,ro_layer,dy_snow,ro_water,tslsnowfall,
c     &  swe_lyr,ro_nsnow,T_old,Tsfc,tsls_threshold,dt,
c     &  melt_flag,soft_snow_d,ro_snow,sum_sprec,sprec_grnd_ml,
c     &  sum_prec,prec,sum_runoff,runoff,snow_d,snow_depth,swe_depth,
c     &  ml_ret)

c Merge layers if the number of layers exceeds some maximum number of
c   layers or if a layer gets thinner than some minimum thickness.
c      CALL MERGE_LAYERS_ML(JJ,ro_layer,dy_snow,swe_lyr,T_old,
c     &  ro_water,max_layers,change_layer,dz_snow_min,melt_flag)
     
c Calculate the temperature of each snow layer.
c      CALL SNOWTEMP_ML(gamma,T_old,Tsfc,JJ,dt,ro_layer,Cp_snow,
c     &  Tf,dy_snow,melt_flag)

c Postprocess the data.
c      CALL POST_PROC_ML(JJ,dy_snow,snow_depth,swe_depth,undef,
c     &  swe_lyr,gamma,ro_layer,melt_flag,T_old,Tf,ro_snow_grid,
c     &  ro_water,xro_snow)

c END J.PFLUG
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GETGAMMA(JJ,ro_layer,gamma)

      implicit none

      include 'snowmodel.inc'

      integer j,JJ
      real ro_layer(nz_max)
      real gamma(nz_max)

c Compute the snow thermal conductivity (gamma) from the snow density.
      do j=1,JJ
        if (ro_layer(j).lt.156.0) then
          gamma(j) = 0.023 + 0.234 * (ro_layer(j)/1000.0)
        else
          gamma(j) = 0.138 - 1.01 * (ro_layer(j)/1000.0) + 3.233 *
     &      (ro_layer(j)/1000.0)**2
        endif
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE POST_PROC_ML(JJ,dy_snow,snow_depth,swe_depth,undef,
     &  swe_lyr,gamma,ro_layer,melt_flag,T_old,Tf,ro_snow_grid,
     &  ro_water,xro_snow)

      implicit none

      include 'snowmodel.inc'

      integer j,JJ

      real dy_snow(nz_max)
      real swe_lyr(nz_max)
      real gamma(nz_max)
      real ro_layer(nz_max)
      real T_old(nz_max)
      integer melt_flag(nz_max)

      real snow_depth,swe_depth,undef,Tf,ro_snow_grid,ro_water,
     &  xro_snow

c Calculate the total snow and swe depth, and the bulk snow density.
      snow_depth = 0.0
      swe_depth = 0.0
      do j=1,JJ
        snow_depth = snow_depth + dy_snow(j)
        swe_depth = swe_depth + swe_lyr(j)
      enddo
      ro_snow_grid = swe_depth * ro_water / snow_depth

c Set any areas outside the snowpack to undef.
      do j=JJ+1,nz_max
        gamma(j) = undef
        ro_layer(j) = undef
        T_old(j) = undef + Tf
        melt_flag(j) = nint(undef)
        dy_snow(j) = undef
        swe_lyr(j) = undef
      enddo

c Clean up the snow density array so there are no values when
c   there is no snow.
      if (snow_depth.eq.0.0) then
        xro_snow = undef
      else
        xro_snow = ro_snow_grid
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MERGE_LAYERS_ML(JJ,ro_layer,dy_snow,swe_lyr,T_old,
     &  ro_water,max_layers,change_layer,dz_snow_min,melt_flag)

      implicit none

      include 'snowmodel.inc'

      integer j,JJ,jjj,j_small,max_layers,icount,k

      real swe_lyr(nz_max)
      real ro_layer(nz_max)
      real dy_snow(nz_max)
      real T_old(nz_max)
      integer melt_flag(nz_max)

      real dy_snow_small,ro_water,dz_snow_min,change_layer

c Merge layers if the number of layers exceeds some maximum number of
c   layers or if a layer gets thinner than some minimum thickness.
c   Do this in snow_depth space because that is the grid the snow
c   temperatures are being calculated on.

c If any layer is thinner than the minimum layer thickness, merge it
c   with the layer below.  If that layer is layer 1, merge it with
c   layer 2.  If there is only one layer left, let it be smaller than
c   the minimum thickness.  Don't do any of this if a new layer is
c   being built; only do it for layers below the top layer.
      change_layer = 0.0

c Count how many layers are less than the minimum thickness, excluding
c   the case where there is only one layer.
      icount = 0
      if (JJ.gt.1) then
c       do j=1,JJ
        do j=1,JJ-1
          if (dy_snow(j).lt.dz_snow_min) then
            icount = icount + 1
          endif
        enddo
      endif



c Note that if two thin layers are together, the merge may take
c   out the other one.
      do k=1,icount
        change_layer = 1.0

c This gets and processes the last occurance.
c       do j=1,JJ
        do j=1,JJ-1
          if (dy_snow(j).lt.dz_snow_min) then
            j_small = j
          endif
        enddo

        if (j_small.eq.1) then
          dy_snow(1) = dy_snow(1) + dy_snow(2)
          swe_lyr(1) = swe_lyr(1) + swe_lyr(2)
          ro_layer(1) = swe_lyr(1) * ro_water / dy_snow(1)
          T_old(1) = T_old(2)
          melt_flag(1) = melt_flag(2)
          JJ = JJ - 1
          do jjj=2,JJ
            dy_snow(jjj) = dy_snow(jjj+1)
              swe_lyr(jjj) = swe_lyr(jjj+1)
            ro_layer(jjj) = swe_lyr(jjj) * ro_water / dy_snow(jjj)
            T_old(jjj) = T_old(jjj+1)
            melt_flag(jjj) = melt_flag(jjj+1)
          enddo
        else
          dy_snow(j_small-1) = dy_snow(j_small-1) + dy_snow(j_small)
          swe_lyr(j_small-1) = swe_lyr(j_small-1) + swe_lyr(j_small)
          ro_layer(j_small-1) = swe_lyr(j_small-1) * ro_water /
     &      dy_snow(j_small-1)
          T_old(j_small-1) = T_old(j_small)
          melt_flag(j_small-1) = melt_flag(j_small)
          JJ = JJ - 1
          do jjj=j_small,JJ
            dy_snow(jjj) = dy_snow(jjj+1)
            swe_lyr(jjj) = swe_lyr(jjj+1)
            ro_layer(jjj) = swe_lyr(jjj) * ro_water / dy_snow(jjj)
            T_old(jjj) = T_old(jjj+1)
            melt_flag(jjj) = melt_flag(jjj+1)
          enddo
        endif
      enddo



c Where the number of layers exceeds some maximum number of layers,
c   find the thinnest layer and merge it with the one below.  For the
c   case where the thinnest layer is the bottom layer, merge it with
c   layer 2.
      if (JJ.eq.max_layers+1) then
        change_layer = 1.0
c Find the thinnest layer.
        dy_snow_small = 1000.0
        do j=1,JJ
          if (dy_snow(j).lt.dy_snow_small) then
            dy_snow_small = dy_snow(j)
            j_small = j
          endif
        enddo

c Adjust accordingly.  Note that layers below the thin layer do not
c   change, unless the thin layer is layer 1.  Also, since the layer
c   is thin, assign the new layer the thick layer temperature.
        if (j_small.eq.1) then
          dy_snow(1) = dy_snow(1) + dy_snow(2)
          swe_lyr(1) = swe_lyr(1) + swe_lyr(2)
          ro_layer(1) = swe_lyr(1) * ro_water / dy_snow(1)
          T_old(1) = T_old(2)
          melt_flag(1) = melt_flag(2)
          JJ = JJ - 1
          do jjj=2,JJ
            dy_snow(jjj) = dy_snow(jjj+1)
              swe_lyr(jjj) = swe_lyr(jjj+1)
            ro_layer(jjj) = swe_lyr(jjj) * ro_water / dy_snow(jjj)
            T_old(jjj) = T_old(jjj+1)
            melt_flag(jjj) = melt_flag(jjj+1)
          enddo
        else
          dy_snow(j_small-1) = dy_snow(j_small-1) + dy_snow(j_small)
          swe_lyr(j_small-1) = swe_lyr(j_small-1) + swe_lyr(j_small)
          ro_layer(j_small-1) = swe_lyr(j_small-1) * ro_water /
     &      dy_snow(j_small-1)
          T_old(j_small-1) = T_old(j_small)
          melt_flag(j_small-1) = melt_flag(j_small)
          JJ = JJ - 1
          do jjj=j_small,JJ
            dy_snow(jjj) = dy_snow(jjj+1)
            swe_lyr(jjj) = swe_lyr(jjj+1)
            ro_layer(jjj) = swe_lyr(jjj) * ro_water / dy_snow(jjj)
            T_old(jjj) = T_old(jjj+1)
            melt_flag(jjj) = melt_flag(jjj+1)
          enddo
        endif
      endif

c Now that we are done with change_layer, set it equal to j_small,
c   the position of the change.
      if (change_layer.eq.1.0) change_layer = real(j_small)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE MELT_SNOW_ML(JJ,swe_lyr,ro_water,ro_layer,Qm,dt,
     &  dy_snow,ro_snowmax,rain,xLf,Cp_snow,Tf,T_old,melt_flag,
     &  runoff,canopy_unload,swe_depth,snow_depth,vegtype,
     &  glacier_melt,sum_glacmelt,sum_swemelt,soft_snow_d,Qe,
     &  sfc_sublim_flag,sum_sfcsublim,swesublim,fc_param,
     &  Cp_water,ml_ret,
     &  change_layer,dz_snow_min,gamma,max_layers,prec,ro_nsnow,
     &  ro_snow,ro_snow_grid,snow_d,sprec_grnd_ml,sum_prec,
     &  sum_runoff,sum_sprec,tsfc,tsls_threshold,tslsnowfall,undef,
     &  xro_snow,Saturn)


      implicit none

      include 'snowmodel.inc'

      integer j,JJ

      real swe_lyr(nz_max)
      real ro_layer(nz_max)
      real dy_snow(nz_max)
      real T_old(nz_max)

      real liqflux,sumup

      real mLayerTheta(nz_max)
      real liqfrac(nz_max)
      real icefrac(nz_max)
      real liqfluxar(nz_max)
      real ml_ret(nz_max)
      real diff(nz_max)
      real Saturn(nz_max)

      real swemelt,extra,ro_water,swe_space,add,runoff,ro_snowmax,
     &  rain,delta_T,xLf,Cp_snow,Tf,dt,Qm,canopy_unload,
     &  potmelt,swe_depth,snow_depth,vegtype,glacier_melt,
     &  sum_glacmelt,sum_swemelt,soft_snow_d,Qe,sfc_sublim_flag,
     &  xLsublim,potsublim,swesublim,sum_sfcsublim,Cp_water

      real fc_param,dt_new
 
      integer melt_flag(nz_max)

      integer max_layers,iteration
      real change_layer,dz_snow_min,prec,ro_nsnow,ro_snow,ro_snow_grid,
     &  snow_d,sprec_grnd_ml,sum_prec,sum_runoff,sum_sprec,tsfc,
     &  tsls_threshold,tslsnowfall,undef,xro_snow
      real gamma(nz_max)

      real swe_depth_old,snow_depth_old,sum_glacmelt_old,runoff_old,
     &  sum_swemelt_old
      real swe_lyr_old(nz_max),dy_snow_old(nz_max),
     &  ro_layer_old(nz_max),melt_flag_old(nz_max)
      integer JJ_old,b
      real sprec_grnd_ml_new,prec_new,runoff_sumup

c Initialize the runoff array.
      runoff = 0.0
      
c SURFACE SUBLIMATION.

c Whether static-surface (non-blowing snow) sublimation is included
c   in the model calculations is controlled by the sfc_sublim_flag.
c   I am waiting for the flux-tower data Matthew and I are collecting
c   in Alaska, to compare with the model simulations, before
c   including this part of the model in all simulations.

c If the sfc_sublim_flag is turned on, the latent heat flux (Qe)
c   calculated in ENBAL is used to add/remove snow from the snowpack.
c   xLsublim = xLf + xLv = 2.5104x10^6 J/kg + 3.334x10^5 J/kg, and
c   potsublim is in m swe.
c      do b = 1,iteration

      if (swe_depth.gt.0.0  .and.  sfc_sublim_flag.eq.1.0) then
        if (Qe.lt.0.0) then
c Compute the snow-surface sublimation (m, swe).
          xLsublim = 2.844e6
          potsublim = (- dt) * Qe / (ro_water * xLsublim)
          swesublim = min(potsublim,swe_depth)
c Save a summing array of the static surface snow sublimation.
          sum_sfcsublim = sum_sfcsublim + swesublim
        else
          swesublim = 0.0
        endif
      else
        swesublim = 0.0
      endif

c Modify the swe layer thicknesses, and reduce the number of layers
c   if needed.
      if (swesublim.gt.0.0) then
c Check to see whether this sublimation requires a layer reduction.
        CALL REDUCE_LAYERS(swesublim,swe_lyr,JJ)

c Build the new snow layer thicknesses, and recalculate the total
c   snow and swe depths.  Assume this sublimated snow does not
c   change the snow density and does not change the soft snow depth.
c   It only reduces the snow depth and the associated swe depth.
        snow_depth = 0.0
        swe_depth = 0.0
        do j=1,JJ
          dy_snow(j) = swe_lyr(j) * ro_water / ro_layer(j)
          snow_depth = snow_depth + dy_snow(j)
          swe_depth = swe_depth + swe_lyr(j)
        enddo
      endif

c J.PFLUG
c Initialize the dynamic timestep iteration
      iteration = 1
      swe_depth_old = swe_depth
      snow_depth_old = snow_depth
      JJ_old = JJ
      swe_lyr_old = swe_lyr
      dy_snow_old = dy_snow
      ro_layer_old = ro_layer
      sum_glacmelt_old = sum_glacmelt
      runoff_old = 0.0
      sum_swemelt_old = sum_swemelt
      melt_flag_old = melt_flag

c if no convergence, return to this point
8888  dt_new = dt/iteration
      runoff_sumup = 0.0

c run through the new number of dynamic timesteps
      do b = 1,iteration
c reinitialize runoff for this dynamic timestep
        runoff = 0.0
c MELTING.

        if (Qm.gt.0.0) then

c Convert the melt energy to water equivalent melt depth (m).
          potmelt = dt_new * Qm / (ro_water * xLf)

c Account for any snowmelt data assimilation.
c       if (icorr_factor_index.lt.0.0) then
c         potmelt_tmp = potmelt * corr_factor
c         swemelt = min(potmelt_tmp,swe_depth)
c Handle the case of no snowmelt data assimilation.
c       else
            swemelt = min(potmelt,swe_depth)
c       endif

c Compute any glacier or permanent snow-field melt (m water equiv.).
          if (vegtype.eq.20.0) then
            glacier_melt = potmelt - swemelt
          else
            glacier_melt = 0.0
          endif

c Save a summing array of the glacier melt.
          sum_glacmelt = sum_glacmelt + glacier_melt

c Save the runoff contribution.
          runoff = runoff + glacier_melt

c Save a summing array of the snow melt.
          sum_swemelt = sum_swemelt + swemelt

c In the presence of melt, zero out the soft snow layer.
          soft_snow_d = 0.0

        else

c These prevent values from the previous time step from being
c   carried through to the next time step.
          swemelt = 0.0
          glacier_melt = 0.0

        endif

c Handle the case where rain and canopy_unload fall on snow-free
c   ground (this is not included in the code below, nor in the
c   PRECIP_ML subroutine, so I include it here).
        if (swe_depth.eq.0.0) then
          runoff = runoff+(rain/iteration)+
     &      (canopy_unload/iteration)
        endif

c Deal with melting snow.

        if (swemelt.gt.0.0) then
c Check to see whether this melt leads to a reduction in layers.
          CALL REDUCE_LAYERS(swemelt,swe_lyr,JJ)

c Build the new snow layer thicknesses, and initiate the melt_flag.
          do j=1,JJ
            dy_snow(j) = swe_lyr(j) * ro_water / ro_layer(j)
            melt_flag(j) = 0
          enddo
        endif

c Add the melt, rain, and canopy unloading (assumed to be wet as rain)
c   to the remaining snow layer thicknesses, up to the maximum snow
c   density, and let the rest of the melt drain out the snowpack bottom
c   as runoff.

c variable to track how much liquid water is available for transport
        if (swe_depth.gt.0.0) then
          sumup = swemelt+(rain/iteration)+
     &      (canopy_unload/iteration)
        
c determine layer properties prior to flux (in fractional form)
          do j=JJ,1,-1
            mLayerTheta(j) = swe_lyr(j)/dy_snow(j)
          enddo

c for each layer, investigate the snowpack composition and then
c instantiate fluxes based on incoming liquid water
          do j=JJ,1,-1

            liqfluxar(j) = sumup

c determine layer ice and liquid content based on initial 
c densities and temperatures
            CALL UPDATE_STATE(T_old(j),mLayerTheta(j),
     &        liqfrac(j),icefrac(j),fc_param,Tf,ro_water,
     &        ml_ret(j))

c determine the amount of liquid water in fractional form
c input in the liquid flux routine
            extra = liqfrac(j)

c determine the unsaturated flux
            CALL SNOWLIQFLX(dt_new,extra,liqfrac(j),icefrac(j),
     &        liqflux,ml_ret(j),Saturn(j))

c check to see if unstable
c unstable is defined as fluxes larger than the available water
            if (liqflux.gt.sumup) then
c if the timestep is larger than 10 minutes and the liquid
c water content is larger than 0.1 mm
              if (dt_new.gt.600.and.sumup.gt.0.0001) then
c increase the timestep and throw message
                iteration = iteration + 1
c                if (iteration.eq.2) then
c                  print *,'employing dynamic timestep for convergence'
c                endif

c reinstantiate state variables
                swe_depth = swe_depth_old
                snow_depth = snow_depth_old
                JJ = JJ_old
                swe_lyr = swe_lyr_old
                dy_snow = dy_snow_old
                ro_layer = ro_layer_old
                sum_glacmelt = sum_glacmelt_old
                runoff = runoff_old
                sum_swemelt = sum_swemelt_old
                melt_flag = melt_flag_old
                goto 8888
              else
c if approaching convergence, force to converge
                liqflux = sumup
              endif
            endif
        
c rebuild the snowpack
            swe_lyr(j) = swe_lyr(j) + sumup - liqflux
            ro_layer(j) = swe_lyr(j) * ro_water/dy_snow(j)
            sumup = liqflux

          enddo

c flux out of the last snow layer goes to runoff
          runoff = runoff + sumup

c recalculate layer depths
          do j=JJ,1,-1
            dy_snow(j) = swe_lyr(j) * ro_water / ro_layer(j)
          enddo

c Also take into account the refreezing of this liquid in a cold
c   snowpack.  Assume that the liquid will fully warm each layer before
c   moving on to the next layer.
          do j=JJ,1,-1

c calculate the amount of liquid that remains in the layer for a step
            if (j.gt.1) then
              diff(j) = liqfluxar(j)-liqfluxar(j-1)
            else
              diff(j) = liqfluxar(j)-runoff
            endif

c Compute the change in temperature that would result if this liquid
c   was used to freeze and raise the snow temperature.
            delta_T = (diff(j) * xLf) / (Cp_snow * dy_snow(j))

c Use this potential temperature change to adjust the snow
c   temperature in the presence of the liquid.
            T_old(j) = T_old(j) + delta_T
            T_old(j) = min(T_old(j),Tf)

c Keep track of which layers have been pushed to Tf.  This will be
c   used in the temperature solution subroutine to fix the snow
c   temperature at Tf (if melt_flag = 1).
            if (T_old(j).eq.Tf) then
              melt_flag(j) = 1
            else
              melt_flag(j) = 0
            endif

          enddo
        endif

c subsample precipitation based on the dynamic timestep
        sprec_grnd_ml_new = sprec_grnd_ml/iteration
        prec_new = prec/iteration

c Account for the accumulation of snow precipitation on the snowpack.
        CALL PRECIP_ML(JJ,ro_layer,dy_snow,ro_water,tslsnowfall,
     &    swe_lyr,ro_nsnow,T_old,Tsfc,tsls_threshold,dt_new,
     &    melt_flag,soft_snow_d,ro_snow,sum_sprec,sprec_grnd_ml_new,
     &    sum_prec,prec_new,sum_runoff,runoff,snow_d,snow_depth,
     &    swe_depth,ml_ret)
c Merge layers if the number of layers exceeds some maximum number of
c   layers or if a layer gets thinner than some minimum thickness.
        CALL MERGE_LAYERS_ML(JJ,ro_layer,dy_snow,swe_lyr,T_old,
     &    ro_water,max_layers,change_layer,dz_snow_min,melt_flag)

c Calculate the temperature of each snow layer.
        CALL SNOWTEMP_ML(gamma,T_old,Tsfc,JJ,dt_new,ro_layer,Cp_snow,
     &    Tf,dy_snow,melt_flag)

c Postprocess the data.
        CALL POST_PROC_ML(JJ,dy_snow,snow_depth,swe_depth,undef,
     &    swe_lyr,gamma,ro_layer,melt_flag,T_old,Tf,ro_snow_grid,
     &    ro_water,xro_snow)


c keep track of the accumulated runoff for the dynamic timestep
      runoff_sumup = runoff_sumup + runoff
      enddo

c total runoff for the set timestep
      runoff = runoff_sumup

c END J.PFLUG
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE REDUCE_LAYERS(swemelt,swe_lyr,JJ)

      implicit none

      include 'snowmodel.inc'

      integer j,JJ
      real swe_lyr(nz_max)
      real eps,swemelt_tmp,swemelt,excess

      eps = 1e-6
      swemelt_tmp = swemelt

c The use of eps here does not allow the vertical grid increment to
c   be less that eps.

      do j=JJ,1,-1
        excess = swe_lyr(j) - swemelt_tmp

c       if (excess.gt.0.0) then
        if (excess.gt.eps) then
          swe_lyr(j) = excess
          JJ = j
          return
c       elseif (excess.eq.0.0) then
        elseif (excess.ge.0.0 .and. excess.le.eps) then
          JJ = j - 1
          return
        else
          swemelt_tmp = - excess
        endif
      enddo

c If there is no snow left.
      JJ = 0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE PRECIP_ML(JJ,ro_layer,dy_snow,ro_water,tslsnowfall,
     &  swe_lyr,ro_nsnow,T_old,Tsfc,tsls_threshold,dt,
     &  melt_flag,soft_snow_d,ro_snow,sum_sprec,sprec_grnd_ml,
     &  sum_prec,prec,sum_runoff,runoff,snow_d,snow_depth,swe_depth,
     &  ml_ret)

      implicit none

      include 'snowmodel.inc'

      integer j,JJ

      real dy_snow(nz_max)
      real ro_layer(nz_max)
      real swe_lyr(nz_max)
      real T_old(nz_max)
      real ml_ret(nz_max)
      real ro_nsnow,ro_water,z_nsnow,tsls_threshold,
     &  z_snowtopl,sweq_topl,Tsfc,tslsnowfall,dt,soft_snow_d,ro_snow,
     &  sum_sprec,sprec_grnd_ml,sum_prec,prec,sum_runoff,runoff,snow_d,
     &  snow_depth,swe_depth

      integer melt_flag(nz_max)

c If the melt from the previous subroutine reduced the snowpack
c   to no snow, reset the time since last snowfall to the threshold,
c   otherwise you will not build a new layer on the bare ground.
      if (JJ.eq.0) tslsnowfall = tsls_threshold

c Create and/or modify the snow c.v.'s to account for new snowfall.
      if (sprec_grnd_ml.gt.0.0) then
        if (tslsnowfall.ge.tsls_threshold) then
c Create a new layer if snowfall has stopped for a period of time
c   greater or equal to the defined threshold.
          JJ = JJ + 1
          z_nsnow = ro_water / ro_nsnow * sprec_grnd_ml
          dy_snow(JJ) = z_nsnow
          ro_layer(JJ) = ro_nsnow
c J.PFLUG
          ml_ret(JJ) = 0.00001
c END J.PFLUG
          swe_lyr(JJ) =  ro_layer(JJ) * dy_snow(JJ) / ro_water
c Define this new snow layer to have the surface temperature.
          T_old(JJ) = Tsfc
          melt_flag(JJ) = 0
        else
c Add to the existing top layer.
          z_nsnow = ro_water / ro_nsnow * sprec_grnd_ml
          z_snowtopl = dy_snow(JJ) + z_nsnow
          sweq_topl = sprec_grnd_ml + dy_snow(JJ) * ro_layer(JJ) /
     &      ro_water
          dy_snow(JJ) = dy_snow(JJ) + z_nsnow
          ro_layer(JJ) = ro_water * sweq_topl / z_snowtopl
          swe_lyr(JJ) = ro_layer(JJ) * dy_snow(JJ) / ro_water
        endif

c Update the total swe and snow depths.
        snow_depth = 0.0
        swe_depth = 0.0
        do j=1,JJ
          snow_depth = snow_depth + dy_snow(j)
          swe_depth = swe_depth + swe_lyr(j)
        enddo
      endif

c Define the time since last snowfall, in hours.  Handle the case
c   where there is no snow on the ground.
      if (sprec_grnd_ml.gt.0.0) then
        tslsnowfall = 0.0
      else
        tslsnowfall = tslsnowfall + dt / 3600.0
      endif
      if (JJ.eq.0) tslsnowfall = tsls_threshold

c The following are set up to be compatible with SnowTran-3D, and
c   are in snow-depth units.  The sum_sprec corrections are done
c   in the SnowTran-3D code.
      soft_snow_d = soft_snow_d + sprec_grnd_ml * ro_water / ro_snow
      snow_d = swe_depth * ro_water / ro_snow
      sum_sprec = sum_sprec + sprec_grnd_ml

c The following are in swe-depth units.
      sum_prec = sum_prec + prec
      sum_runoff = sum_runoff + runoff

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE DDENSITY_ML(ro_layer,Tf,dt,ro_water,ro_ice,
     &  T_old,JJ,dy_snow,A1,A2)

      implicit none

      include 'snowmodel.inc'

      integer j,JJ,jjj

      real dy_snow(nz_max)
      real ro_layer(nz_max)
      real T_old(nz_max)
      real sweqstar(nz_max)
      real sweql(nz_max)

      real A1,A2,ro_water,ro_ice,dt,Tf

      if (JJ.gt.0) then

        do j=1,JJ
          sweql(j) = ro_layer(j) / ro_water * dy_snow(j)
        enddo

        do jjj=1,JJ
          sweqstar(jjj) = sweql(jjj) / 2.0
          do j=jjj+1,JJ
            sweqstar(jjj) = sweqstar(jjj) + sweql(j)
          enddo
        enddo

        do j=1,JJ
          ro_layer(j) = ro_layer(j) + dt * (A1 * sweqstar(j) *
     &      ro_layer(j) *
     &      exp(-0.08*(Tf-T_old(j))) * exp(-A2*ro_layer(j)))
          ro_layer(j) = min(ro_ice,ro_layer(j))
          dy_snow(j) = sweql(j) * ro_water / ro_layer(j)
        enddo

      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SNOWTEMP_ML(gamma,T_old,Tsfc,JJ,dt,ro_layer,Cp_snow,
     &  Tf,dy_snow,melt_flag)

      implicit none

      include 'snowmodel.inc'

      real gamma(nz_max)
      real ro_layer(nz_max)
      real dy_snow(nz_max)
      real g_b_ns(nz_max+1)
      real f_n(nz_max+1)
      real aN(nz_max)
      real aP0(nz_max)
      real aS(nz_max)
      real dely_p(nz_max+1)
      real dy_p(nz_max)
      real y_crds(nz_max+2)
      real y_wall(nz_max+1)
      real A_sub(nz_max)
      real A_super(nz_max)
      real A_main(nz_max)
      real b_vector(nz_max)
      real T_old(nz_max)
      real Sc(nz_max)
      real Sp(nz_max)

      integer melt_flag(nz_max)

      integer j,JJ
      real Tsfc,T_N,bc_N,bc_S,Cp_snow,Tf,dt,Tsg

c Define the snow thermal conductivity (gamma) for each layer.
      CALL GETGAMMA(JJ,ro_layer,gamma)

      if (JJ.gt.1) then

c Update the control volume information.
        CALL GETCV(JJ,dy_p,dy_snow)
        CALL CV_INFO(dely_p,f_n,y_crds,y_wall,dy_p,JJ)

c Compute the general equation coefficients.
        CALL GAMMA1(g_b_ns,gamma,f_n,JJ)
        CALL GE_COEF(aN,aS,aP0,dy_p,dely_p,g_b_ns,dt,JJ,
     &    ro_layer,Cp_snow)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
c Account for the boundary conditions.
c   South boundary condition:
c     For T_S = known, define 
c       bc_S = aS(1) * T_S;         where T_S = known
c     For dT_S/dn = 0, define
c       bc_S = 0.0
c       aS(1) = 0.0
c   North boundary condition:
c     For T_N = known, define 
c       bc_N = aN(JJ) * T_N;        where T_N = known
c     For dT_N/dn = 0, define
c       bc_N = 0.0
c       aN(JJ) = 0.0
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c Define the upper and lower boundary conditions.
        T_N = Tsfc
        bc_N = aN(JJ) * T_N
        bc_S = 0.0
        aS(1) = 0.0

c Provide the source terms.

c Force the source terms to produce Tf at the positions where melting
c   occurred during this time step.
        do j=1,JJ
          if (melt_flag(j).eq.1) then
            Sc(j) = 10e30 * Tf
            Sp(j) = -10e30
          else
            Sc(j) = 0.0
            Sp(j) = 0.0
          endif
        enddo

c Configure the information for the matrix solver.
        CALL PREPSOLVE(A_sub,A_super,A_main,b_vector,T_old,
     &    dy_p,bc_S,bc_N,Sc,Sp,aN,aS,aP0,JJ)

c Solve the system of equations.
        CALL TRISOLVE(T_old,A_sub,A_main,A_super,b_vector,JJ)

      elseif (JJ.eq.1) then
c Assume that the snow-ground interface temperature is -1.0 C.
c J.PFLUG
c changed for melt-out conditions
        if (melt_flag(1).eq.1) then
          T_old(1) = Tf
        else
          Tsg = Tf - 1.0
          T_old(1) = 0.5 * (Tsg + Tsfc)
        endif
c END J.PFLUG
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GETCV(JJ,dy_p,dy_snow)

      implicit none

      include 'snowmodel.inc'

      real dy_p(nz_max)
      real dy_snow(nz_max)

      integer j,JJ

c Provide values of Control Volume size in the y direction.
      do j=1,JJ
        dy_p(j) = dy_snow(j)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE PREPSOLVE(A_sub,A_super,A_main,b_vector,T_old,
     &  dy_p,bc_S,bc_N,Sc,Sp,aN,aS,aP0,JJ)

      implicit none

      include 'snowmodel.inc'

      real aP(nz_max)
      real aN(nz_max)
      real aS(nz_max)
      real Sp(nz_max)
      real Sc(nz_max)
      real aP0(nz_max)
      real dy_p(nz_max)
      real T_old(nz_max)
      real b_vector(nz_max)
      real A_sub(nz_max)
      real A_super(nz_max)
      real A_main(nz_max)

      integer j,jj
      real bc_S,bc_N

c Compute matrix diagonal and b coeffs.
      do j=1,JJ
        aP(j) = aN(j) + aS(j) + aP0(j) - Sp(j) * dy_p(j)
        b_vector(j) = Sc(j) * dy_p(j) + aP0(j) * T_old(j)
      enddo

c Modify b to account for dirichlet boundary conditions.
      b_vector(1) = b_vector(1) + bc_S
      b_vector(JJ) = b_vector(JJ) + bc_N

c Prepare to call the tridiagonal solver.
      do j=1,JJ-1
        A_sub(j) = - aS(j+1)
        A_super(j) = - aN(j)
      enddo

      do j=1,JJ
        A_main(j) = aP(j)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CV_INFO(dely_p,f_n,y_crds,y_wall,dy_p,JJ)

      implicit none

      include 'snowmodel.inc'

      real dy_pbc(nz_max+2)
      real dely_p(nz_max+1)
      real f_n(nz_max+1)
      real dy_p(nz_max)
      real y_crds(nz_max+2)
      real y_wall(nz_max+1)

      integer j,JJ
      real temp

c PRESSURE CONTROL VOLUME SIZE AND POSITION INFORMATION

c Include exterior boundary pressure grid points.
      dy_pbc(1) = 0.0
      do j=2,JJ+1
        dy_pbc(j) = dy_p(j-1)
      enddo
      dy_pbc(JJ+2) = 0.0

c Compute the distance between pressure grid points.
      do j=1,JJ+1
        dely_p(j) = .5 * (dy_pbc(j) + dy_pbc(j+1))
      enddo

c Compute the distance between the pressure grid points and the control
c   volume wall.  (The following is true because the grid points do
c   pressure are defined to be in the center of the control volume.)
c   And then compute f_e and f_n.  These two steps are combined below.
      do j=1,JJ+1
        f_n(j) = .5 * dy_pbc(j+1) / dely_p(j)
      enddo

c Compute the x and y coordinates of the pressure c.v. grid points,
c   including boundaries.
      temp = 0.0
      do j=1,JJ+2
        y_crds(j) = temp + .5 * dy_pbc(j)
        temp = temp + dy_pbc(j)
      enddo

c Compute the x and y coordinates of the pressure c.v. walls.
      y_wall(1) = 0.0
      do j=2,JJ+1
        y_wall(j) = y_wall(j-1) + dy_p(j-1)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GAMMA1(g_b_ns,gamma,f_n,JJ)

      implicit none

      include 'snowmodel.inc'

      real g_b_ns(nz_max+1)
      real gamma(nz_max)
      real g_ns(nz_max+2)
      real f_n(nz_max+1)

      integer j,JJ

c This provides gamma information on c.v. walls.

c Include gamma just outside of n, s boundaries.
      g_ns(1) = gamma(1)
      do j=2,JJ+1
        g_ns(j) = gamma(j-1)
      enddo
      g_ns(JJ+2) = gamma(JJ)

c Compute gamma (diffusion coefficient) at the n, s control
c   volume boundaries using equation 4.9, p. 45.
      do j=1,JJ+1
        g_b_ns(j) = 1.0/((1.0 - f_n(j))/g_ns(j) + f_n(j)/g_ns(j+1))
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GE_COEF(aN,aS,aP0,dy_p,dely_p,g_b_ns,dt,JJ,
     &  ro_layer,Cp_snow)

      implicit none

      include 'snowmodel.inc'

      real aN(nz_max)
      real aS(nz_max)
      real aP0(nz_max)
      real dely_p(nz_max+1)
      real g_b_ns(nz_max+1)
      real dy_p(nz_max)
      real ro_layer(nz_max)

      integer j,JJ
      real Cp_snow,dt

c CALCULATE THE COEFFICIENTS aP, for the general phi equation.
      do j=2,JJ+1
        aN(j-1) = g_b_ns(j)   / dely_p(j)
        aS(j-1) = g_b_ns(j-1) / dely_p(j-1)
      enddo

      do j=1,JJ
        aP0(j) = ro_layer(j) * Cp_snow * dy_p(j) / dt
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE TRISOLVE(x,asub,amain,asuper,b,JJ)

      implicit none

      include 'snowmodel.inc'

      real asub(nz_max)
      real asuper(nz_max)
      real amain(nz_max)
      real b(nz_max)
      real x(nz_max)
      real z(nz_max)
      real lmain(nz_max)
      real lsub(nz_max)
      real usuper(nz_max)

      integer j,JJ

      lmain(1) = amain(1)
      usuper(1) = asuper(1)/lmain(1)

      do j=2,JJ-1
        lsub(j-1) = asub(j-1)
        lmain(j) = amain(j) - lsub(j-1) * usuper(j-1)
        usuper(j) = asuper(j) / lmain(j)
      enddo

      lsub(JJ-1) = asub(JJ-1)
      lmain(JJ) = amain(JJ) - lsub(JJ-1) * usuper(JJ-1)
      z(1) = b(1) / lmain(1)

      do j=2,JJ
        z(j) = 1.0 / lmain(j) * (b(j) - lsub(j-1) * z(j-1))
      enddo

      x(JJ) = z(JJ)

      do j=JJ-1,1,-1
        x(j) = z(j) - usuper(j) * x(j+1)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ZERO_SNOW(nx,ny,snow_depth,ro_snow_grid,ro_snow,
     &  swe_depth,swe_depth_old,canopy_int_old,JJ,sum_swemelt,
     &  tslsnowfall,dy_snow,swe_lyr,ro_layer,T_old,sum_sprec,
     &  multilayer_snowpack,tsls_threshold)

      implicit none

      include 'snowmodel.inc'

      integer nx,ny,i,j,k

      integer multilayer_snowpack
      integer JJ(nx_max,ny_max)

      real tsls_threshold,ro_snow
      real tslsnowfall(nx_max,ny_max)
      real dy_snow(nx_max,ny_max,nz_max)
      real swe_lyr(nx_max,ny_max,nz_max)
      real ro_layer(nx_max,ny_max,nz_max)
      real T_old(nx_max,ny_max,nz_max)
      real swe_depth_old(nx_max,ny_max)
      real canopy_int_old(nx_max,ny_max)
      real swe_depth(nx_max,ny_max)
      real snow_depth(nx_max,ny_max)
      real ro_snow_grid(nx_max,ny_max)
      real sum_sprec(nx_max,ny_max)
      real sum_swemelt(nx_max,ny_max)

      print *,'ZEROING OUT THE SNOW ARRAYS'
      print *,'ZEROING OUT THE SNOW ARRAYS'
      print *,'ZEROING OUT THE SNOW ARRAYS'

      do j=1,ny
        do i=1,nx
          canopy_int_old(i,j) = 0.0
          swe_depth_old(i,j) = 0.0
          snow_depth(i,j) = 0.0
          ro_snow_grid(i,j) = ro_snow
          swe_depth(i,j) = 0.0
          sum_sprec(i,j) = 0.0
          sum_swemelt(i,j) = 0.0
        enddo
      enddo

      if (multilayer_snowpack.eq.1) then
        do j=1,ny
          do i=1,nx
            tslsnowfall(i,j) = tsls_threshold
            do k=1,JJ(i,j)
              dy_snow(i,j,k) = 0.0
              swe_lyr(i,j,k) = 0.0
              ro_layer(i,j,k) = ro_snow
              T_old(i,j,k) = 273.16
            enddo
            JJ(i,j) = 0
          enddo
        enddo
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ZERO_SEAICE_SNOW(nx,ny,snow_depth,ro_snow_grid,
     &  ro_snow,swe_depth,swe_depth_old,canopy_int_old,JJ,
     &  tslsnowfall,dy_snow,swe_lyr,ro_layer,T_old,
     &  multilayer_snowpack,tsls_threshold,seaice_conc)

      implicit none

      include 'snowmodel.inc'

      integer nx,ny,i,j,k

      integer multilayer_snowpack
      integer JJ(nx_max,ny_max)

      real tsls_threshold,ro_snow
      real tslsnowfall(nx_max,ny_max)
      real dy_snow(nx_max,ny_max,nz_max)
      real swe_lyr(nx_max,ny_max,nz_max)
      real ro_layer(nx_max,ny_max,nz_max)
      real T_old(nx_max,ny_max,nz_max)
      real swe_depth_old(nx_max,ny_max)
      real canopy_int_old(nx_max,ny_max)
      real swe_depth(nx_max,ny_max)
      real snow_depth(nx_max,ny_max)
      real ro_snow_grid(nx_max,ny_max)
      real seaice_conc(nx_max,ny_max)

      do j=1,ny
        do i=1,nx
          if (seaice_conc(i,j).eq.0.0) then
            canopy_int_old(i,j) = 0.0
            swe_depth_old(i,j) = 0.0
            snow_depth(i,j) = 0.0
            ro_snow_grid(i,j) = ro_snow
            swe_depth(i,j) = 0.0
          endif
        enddo
      enddo

      if (multilayer_snowpack.eq.1) then
        do j=1,ny
          do i=1,nx
            if (seaice_conc(i,j).eq.0.0) then
              tslsnowfall(i,j) = tsls_threshold
              do k=1,JJ(i,j)
                dy_snow(i,j,k) = 0.0
                swe_lyr(i,j,k) = 0.0
                ro_layer(i,j,k) = ro_snow
                T_old(i,j,k) = 273.16
              enddo
              JJ(i,j) = 0
            endif
          enddo
        enddo
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SNOWLIQFLX(dt,mLayerVolFracLiqTrial,
     &  liqfrac,icefrac,liqflux,
     &  mLayerThetaResid,relSaturn)

      implicit none

c JPFLUG
c public variables

      real mLayerVolFracLiqTrial,liqfrac,icefrac,liqflux,dt

      real pore_space,availCap,relSaturn,iden_ice,mLayerThetaResid,
     &  k_snow,mult,mw

      parameter (iden_ice = 917)

c compute the pore space. Note that water can fill some of this
c pore space
        pore_space = 1.0 - icefrac

c compute fluxes
c check that flow occurs
        if(mLayerVolFracLiqTrial>mLayerThetaResid)then

c compute the available capacity
          availCap = pore_space - mLayerThetaResid

c compute the saturated hydraulic conductivity
c method adapted from Colbeck (1972).
          k_snow = (3.41875e-5)*exp(15.9*(1-icefrac))
         
c compute the relative saturation
          if (availCap.gt.0.0) then
            relSaturn = (mLayerVolFracLiqTrial -
     &        mLayerThetaResid) / availCap
          else
            relSaturn = 1.0
          endif

c calculate the flux out of the snowpack or layer
           liqflux = dt*k_snow*relSaturn**3.0
        else
            
c flow does not occur
          liqflux = 0.0
          relSaturn = 0.0

        endif

c END JPFLUG
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE UPDATE_STATE(Tk,mLayerTheta,mLayerVolFracLiq,
     &  mLayerVolFracIce,fc_param,Tf,ro_water,
     &  retent)

      implicit none

c JPFLUG
      real fracliquid,fc_param,Tf,Tk,mLayerVolFracLiq,
     &  mLayerVolFracIce,iden_ice,mLayerTheta,
     &  ro_water,retent

      parameter (iden_ice = 917)

c fraction of liquid water
      fracliquid = 1.0/(1.0 + (fc_param*(Tf -
     &  min(Tk,Tf)))**2.0)

c calculate volumetric fractions of liquid and ice content
      mLayerVolFracLiq = fracliquid*mLayerTheta
      mLayerVolFracIce = (1.0 - fracliquid)*mLayerTheta*
     &  (ro_water/iden_ice)

c calculate the retention
      retent = min(0.02,max(retent,0.75*mLayerVolFracLiq))

c END JPFLUG
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
