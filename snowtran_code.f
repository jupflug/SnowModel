c snowtran_code.f

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc     Snow-Transport Modeling System - 3D (SnowTran-3D)    cccccc
ccccc                    Copyright (C) 1998                    cccccc
ccccc          by Glen E. Liston, InterWorks Consulting        cccccc
ccccc                    All Rights Reserved                   cccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c This FORTRAN code receives inputs of wind speed, wind direction,
c   air temperature, relative humidity, vegetation type, topography,
c   and precipitation, and it outputs snow depth, saltation flux,
c   suspended flux, sublimation of blowing snow, and the snow depth
c   changes resulting from these processes.
c
c All units are in m, kg, s, K.
c
c This model is described in the paper:
c   A Snow-Transport Model for Complex Terrain, by Glen E. Liston
c   and Matthew Sturm, Journal of Glaciology, 1998, Vol. 44,
c   No. 148, pages 498-516.
c
c The author of this code is:
c   Dr. Glen E. Liston
c   InterWorks Consulting
c   15621 SnowMan Road
c   Loveland, Colorado 80538
c
c To run in 2-D mode, set nx = 3 and look at the data at i = 2.
c   This is required because of the boundary conditions imposed
c   along i = 1 and 3.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine SNOWTRAN_CODE(bc_flag,bs_flag,C_z,
     &  conc_salt,deltax,deltay,dh_salt,dh_salt_u,dh_salt_v,
     &  dh_susp,dh_susp_u,dh_susp_v,dt,dz_susp,fall_vel,fetch,
     &  gravity,h_const,h_star,ht_rhobs,ht_windobs,index_ue,
     &  index_uw,index_vn,index_vs,iter,nx,ny,pi,Qsalt,Qsalt_max,
     &  Qsalt_maxu,Qsalt_maxv,Qsalt_u,Qsalt_v,Qsubl,Qsusp,
     &  Qsusp_u,Qsusp_v,rh_grid,ro_air,ro_snow,ro_water,snow_d,
     &  snow_d_init,snow_z0,soft_snow_d,sprec,sum_glacmelt,
     &  subgrid_flag,wbal_salt,wbal_susp,wbal_qsubl,sum_sprec,
     &  tabler_ee,tabler_ne,tabler_nn,tabler_nw,tabler_se,
     &  tabler_ss,tabler_sw,tabler_ww,tair_grid,topo,topo_land,
     &  topoflag,tp_scale,twolayer_flag,Up_const,Ur_const,Utau,
     &  Utau_t,uwind_grid,veg_z0,vegsnowd_xy,vegtype,vonKarman,
     &  vwind_grid,wind_min,winddir_flag,winddir_grid,
     &  windspd_flag,windspd_grid,xmu,z_0,ztop_susp,erosion_dist,
     &  run_enbal,run_snowpack,wbal_subgrid,sum_qsubl,sum_trans,
     &  swe_depth,snow_depth,ro_snow_grid,sum_prec,sum_runoff,
     &  sum_Qcs,canopy_int,w_balance,sum_sfcsublim,tabler_dir,
     &  slope_adjust,Utau_t_const,Utau_t_flag,ro_soft_snow_old,
     &  ro_soft_snow,ro_nsnow,prec,Qcs,runoff,d_canopy_int,
     &  glacier_melt,swe_depth_old,swesublim,canopy_unload,
     &  canopy_int_old,iter_start,multilayer_snowpack,swe_lyr,
     &  JJ,dy_snow,ro_layer,curve_lg_scale_flag,curve_wt_lg,
     &  seaice_run,seaice_conc,tslsnowfall,T_old,tsls_threshold)

      implicit none

      include 'snowmodel.inc'

      integer iter,nx,ny,i,j,iter_start

      real ro_snow,ro_water,ro_air,gravity,vonKarman,snow_z0
      real deltax,deltay,dt
      real fetch,xmu,C_z,h_const,wind_min,windspd_flag
      real Up_const,dz_susp,ztop_susp,fall_vel,Ur_const
      real Utau_t_const,pi,bc_flag,topoflag,Utau_t_flag
      real ht_windobs,ht_rhobs,bs_flag,twolayer_flag
      real subgrid_flag,tp_scale,winddir_flag,erosion_dist
      real run_enbal,run_snowpack,tabler_dir,slope_adjust

      real topo_land(nx_max,ny_max)
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

      real tabler_nn_orig(nx_max,ny_max)
      real tabler_ss_orig(nx_max,ny_max)
      real tabler_ee_orig(nx_max,ny_max)
      real tabler_ww_orig(nx_max,ny_max)
      real tabler_ne_orig(nx_max,ny_max)
      real tabler_se_orig(nx_max,ny_max)
      real tabler_sw_orig(nx_max,ny_max)
      real tabler_nw_orig(nx_max,ny_max)
      real snow_d_tabler(nx_max,ny_max)
      real topo_tmp(nx_max,ny_max)

      real uwind_grid(nx_max,ny_max),vwind_grid(nx_max,ny_max)
      real windspd_grid(nx_max,ny_max),winddir_grid(nx_max,ny_max)
      real tair_grid(nx_max,ny_max),sprec(nx_max,ny_max)
      real rh_grid(nx_max,ny_max)

      integer index_ue(ny_max,2*nx_max+1),index_uw(ny_max,2*nx_max+1)
      integer index_vn(nx_max,2*ny_max+1),index_vs(nx_max,2*ny_max+1)

      real snow_d(nx_max,ny_max)
      real snow_depth(nx_max,ny_max)
      real swe_depth(nx_max,ny_max)
      real ro_snow_grid(nx_max,ny_max)
      real ro_soft_snow(nx_max,ny_max)
      real ro_soft_snow_old(nx_max,ny_max)
      real ro_nsnow(nx_max,ny_max)
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

      real dh_subgrid(nx_max,ny_max)
      real Qsubl(nx_max,ny_max)

      real sum_sprec(nx_max,ny_max)
      real wbal_qsubl(nx_max,ny_max)
      real wbal_salt(nx_max,ny_max)
      real wbal_susp(nx_max,ny_max)
      real wbal_subgrid(nx_max,ny_max)
      real sum_qsubl(nx_max,ny_max)
      real sum_trans(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)

      real prec(nx_max,ny_max)
      real Qcs(nx_max,ny_max)
      real runoff(nx_max,ny_max)
      real d_canopy_int(nx_max,ny_max)
      real glacier_melt(nx_max,ny_max)
      real swe_depth_old(nx_max,ny_max)
      real swesublim(nx_max,ny_max)
      real canopy_unload(nx_max,ny_max)
      real canopy_int_old(nx_max,ny_max)

      real vegsnowd_xy(nx_max,ny_max)
      real veg_z0(nx_max,ny_max)

      real sum_glacmelt(nx_max,ny_max),w_balance(nx_max,ny_max),
     &  sum_prec(nx_max,ny_max),sum_runoff(nx_max,ny_max),
     &  sum_Qcs(nx_max,ny_max),canopy_int(nx_max,ny_max),
     &  sum_sfcsublim(nx_max,ny_max)

      integer multilayer_snowpack,k
      integer JJ(nx_max,ny_max)
      real swe_change_tmp,swe_change,tsls_threshold
      real swe_lyr_z(nz_max)
      real swe_lyr(nx_max,ny_max,nz_max)
      real dy_snow(nx_max,ny_max,nz_max)
      real ro_layer(nx_max,ny_max,nz_max)
      real T_old(nx_max,ny_max,nz_max)
      real tslsnowfall(nx_max,ny_max)

      real curve_lg_scale_flag
      real curve_wt_lg(nx_max,ny_max)

      real seaice_run
      real seaice_conc(nx_max,ny_max)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Perform some intialization steps that are unique to SnowTran-3D.
      if (iter.eq.iter_start) then

        print *,
     & 'cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
        print *,
     & 'c     Snow-Transport Modeling System - 3D (SnowTran-3D)    c'
        print *,
     & 'c                    Copyright (C) 1998                    c'
        print *,
     & 'c          by Glen E. Liston, InterWorks Consulting        c'
        print *,
     & 'c                    All Rights Reserved                   c'
        print *,
     & 'cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'

        if (subgrid_flag.eq.1.0) then

c Check to make sure topoflag = 0.0.
          if (topoflag.eq.1.0) then
            print *,'If subgrid_flag=1.0, then topoflag must = 0.0'
            print *,'  Correct this in snowmodel.par to continue.'
            stop
          endif

c The Tabler surfaces were originally developed assuming the grid
c   increment would never be less than 1.0 m.  You probably should
c   not run it with deltax and deltay less than 1.0 without some
c   further testing.
          if (deltax.lt.1.0 .or. deltay.lt.1.0) then
            print *,'The Tabler subgrid algorithm has not been'
            print *,'tested for deltax and/or deltay less than'
            print *,'1.0 m.  You should probably do some testing'
            print *,'before running the model at less than 1.0-m'
            print *,'resolution.  Acually I am pretty sure it will'
            print *,'run, but it will not generate the correct'
            print *,'snow-depth profiles.'
            stop
          endif

c If this is the first time through, generate the Tabler snow
c   accumulation surfaces for the land topography.
          call tabler_3d(nx,ny,topo_land,deltax,deltay,
     &      tabler_ww_orig,tabler_ee_orig,tabler_ss_orig,
     &      tabler_nn_orig,erosion_dist,tabler_ne_orig,
     &      tabler_se_orig,tabler_sw_orig,tabler_nw_orig,
     &      slope_adjust)

c As part of generating the Tabler surfaces, a -8888.0 has been
c   used to identify the areas immediately upwind of any Tabler
c   drift trap that is an erosion area where no snow is allowed
c   to accumulate (see 'erosion_dist' in snowmodel.par).  Now
c   take those areas and set them equal to the snow-holding depth
c   to keep the snow relatively thin in those areas.
          do i=1,nx
            do j=1,ny
              tabler_nn_orig(i,j) =
     &          max(tabler_nn_orig(i,j),vegsnowd_xy(i,j))
              tabler_ne_orig(i,j) =
     &          max(tabler_ne_orig(i,j),vegsnowd_xy(i,j))
              tabler_ee_orig(i,j) =
     &          max(tabler_ee_orig(i,j),vegsnowd_xy(i,j))
              tabler_se_orig(i,j) =
     &          max(tabler_se_orig(i,j),vegsnowd_xy(i,j))
              tabler_ss_orig(i,j) =
     &          max(tabler_ss_orig(i,j),vegsnowd_xy(i,j))
              tabler_sw_orig(i,j) =
     &          max(tabler_sw_orig(i,j),vegsnowd_xy(i,j))
              tabler_ww_orig(i,j) =
     &          max(tabler_ww_orig(i,j),vegsnowd_xy(i,j))
              tabler_nw_orig(i,j) =
     &          max(tabler_nw_orig(i,j),vegsnowd_xy(i,j))
            enddo
          enddo

c If you don't want to write out these distributions, comment out
c   the following lines.
          open(51,file='outputs/tabler_sfcs.gdat',
     &      form='unformatted',access='direct',recl=4*nx*ny)
          write(51,rec=1) ((tabler_ww_orig(i,j),i=1,nx),j=1,ny)
          write(51,rec=2) ((tabler_ee_orig(i,j),i=1,nx),j=1,ny)
          write(51,rec=3) ((tabler_ss_orig(i,j),i=1,nx),j=1,ny)
          write(51,rec=4) ((tabler_nn_orig(i,j),i=1,nx),j=1,ny)
          write(51,rec=5) ((tabler_ne_orig(i,j),i=1,nx),j=1,ny)
          write(51,rec=6) ((tabler_se_orig(i,j),i=1,nx),j=1,ny)
          write(51,rec=7) ((tabler_sw_orig(i,j),i=1,nx),j=1,ny)
          write(51,rec=8) ((tabler_nw_orig(i,j),i=1,nx),j=1,ny)
          write(51,rec=9) ((topo_land(i,j),i=1,nx),j=1,ny)
          close (51)

c This available if you want to save a specific Tabler surface at
c   each time step.
c         open(52,file='outputs/tabler_sfcs_iter.gdat',
c    &      form='unformatted',access='direct',recl=4*nx*ny)

        endif

      endif

c Print out some basic run information to the screen.
      print 102, windspd_flag,winddir_flag
  102 format(25x,'    wind spd = ',f5.2,'   wind dir = ',f4.0)

      if (subgrid_flag.eq.1.0) then

c Generate the tabler surfaces at this time step, assuming the
c   snow surface is the topographic surface.

c Take the snow depth coming out of SnowPack, and before any wind
c   redistribution has been applied at this time step, and add it
c   to topo_land.  Then use this to create the Tabler surface that
c   will be used for this time step.  Also define this depth to be
c   dependent on the SnowPack spatially distributed snow density,
c   not the constant density used in SnowTran.
        do i=1,nx
          do j=1,ny
            snow_d_tabler(i,j) = swe_depth(i,j) *
     &        ro_water / ro_snow_grid(i,j)
            topo_tmp(i,j) = tp_scale * snow_d_tabler(i,j) +
     &        topo_land(i,j)
          enddo
        enddo

        call tabler_3d(nx,ny,topo_tmp,deltax,deltay,
     &    tabler_ww_orig,tabler_ee_orig,tabler_ss_orig,
     &    tabler_nn_orig,erosion_dist,tabler_ne_orig,
     &    tabler_se_orig,tabler_sw_orig,tabler_nw_orig,
     &    slope_adjust)

c The Tabler surfaces that were just generated have had topo_tmp
c   subtracted off of them, giving just the drift profiles with
c   things like zero drift depth on ridges and windwards slopes.
c   So, add the snow depth, prior to any wind redistribution, to
c   these Tabler surfaces.  This will be the maximum snow depth
c   allowed as part of the wind redistribution.
        do i=1,nx
          do j=1,ny
            tabler_ww(i,j) = tp_scale * snow_d_tabler(i,j) +
     &        tabler_ww_orig(i,j)
            tabler_ee(i,j) = tp_scale * snow_d_tabler(i,j) +
     &        tabler_ee_orig(i,j)
            tabler_ss(i,j) = tp_scale * snow_d_tabler(i,j) +
     &        tabler_ss_orig(i,j)
            tabler_nn(i,j) = tp_scale * snow_d_tabler(i,j) +
     &        tabler_nn_orig(i,j)
            tabler_ne(i,j) = tp_scale * snow_d_tabler(i,j) +
     &        tabler_ne_orig(i,j)
            tabler_se(i,j) = tp_scale * snow_d_tabler(i,j) +
     &        tabler_se_orig(i,j)
            tabler_sw(i,j) = tp_scale * snow_d_tabler(i,j) +
     &        tabler_sw_orig(i,j)
            tabler_nw(i,j) = tp_scale * snow_d_tabler(i,j) +
     &        tabler_nw_orig(i,j)
          enddo
        enddo

c As part of generating the Tabler surfaces, a -8888.0 has been
c   used to identify the areas immediately upwind of any Tabler
c   drift trap that is an erosion area where no snow is allowed
c   to accumulate (see 'erosion_dist' in snowmodel.par).  Now
c   take those areas and set them equal to the snow-holding depth
c   to keep the snow relatively thin in those areas.
        do i=1,nx
          do j=1,ny
            tabler_nn(i,j) = max(tabler_nn(i,j),vegsnowd_xy(i,j))
            tabler_ne(i,j) = max(tabler_ne(i,j),vegsnowd_xy(i,j))
            tabler_ee(i,j) = max(tabler_ee(i,j),vegsnowd_xy(i,j))
            tabler_se(i,j) = max(tabler_se(i,j),vegsnowd_xy(i,j))
            tabler_ss(i,j) = max(tabler_ss(i,j),vegsnowd_xy(i,j))
            tabler_sw(i,j) = max(tabler_sw(i,j),vegsnowd_xy(i,j))
            tabler_ww(i,j) = max(tabler_ww(i,j),vegsnowd_xy(i,j))
            tabler_nw(i,j) = max(tabler_nw(i,j),vegsnowd_xy(i,j))
          enddo
        enddo

c The following saves the calculated Tabler surface at each time
c   step.  You can comment this out if you don't want to write
c   them out.
c       if (tabler_dir.eq.0.0) then
c         write(52,rec=iter) ((tabler_nn(i,j),i=1,nx),j=1,ny)
c       elseif (tabler_dir.eq.45.0) then
c         write(52,rec=iter) ((tabler_ne(i,j),i=1,nx),j=1,ny)
c       elseif (tabler_dir.eq.90.0) then
c         write(52,rec=iter) ((tabler_ee(i,j),i=1,nx),j=1,ny)
c       elseif (tabler_dir.eq.135.0) then
c         write(52,rec=iter) ((tabler_se(i,j),i=1,nx),j=1,ny)
c       elseif (tabler_dir.eq.180.0) then
c         write(52,rec=iter) ((tabler_ss(i,j),i=1,nx),j=1,ny)
c       elseif (tabler_dir.eq.225.0) then
c         write(52,rec=iter) ((tabler_sw(i,j),i=1,nx),j=1,ny)
c       elseif (tabler_dir.eq.270.0) then
c         write(52,rec=iter) ((tabler_ww(i,j),i=1,nx),j=1,ny)
c       elseif (tabler_dir.eq.315.0) then
c         write(52,rec=iter) ((tabler_nw(i,j),i=1,nx),j=1,ny)
c       endif

      endif

c In SnowTran-3D, the summed snow precipitation must be in units
c   of snow-depth.  The rest of the routines assume that it is in
c   swe units.
      do j=1,ny
        do i=1,nx
          sum_sprec(i,j) = sum_sprec(i,j) * ro_water / ro_snow
        enddo
      enddo

c If running EnBal and SnowPack, then don't need to run these
c   two routines. 
      if (run_enbal.ne.1.0 .and. run_snowpack.ne.1.0) then
        print *,'I am not sure you can configure the model like'
        print *,'this anymore.  You may need to always run'
        print *,'SnowTran-3D with SnowPack now.'
        stop
c Add the new precipitation to the snowpack.
        call precip(snow_d,sprec,nx,ny,ro_snow,
     &    ro_water,sum_sprec,soft_snow_d,sum_prec)

c If the two-layer scheme is turned on, update the thicknesses
c   of the hard and soft layers.
        if (twolayer_flag.eq.1.0) then
          call twolayer1(nx,ny,soft_snow_d,tair_grid)
        endif
      endif

c Update the threshold friction velocity.
      if (Utau_t_flag.eq.0.0) then
        if (curve_lg_scale_flag.eq.1.0) then
          do j=1,ny
            do i=1,nx
              Utau_t(i,j) = curve_wt_lg(i,j) * Utau_t_const
            enddo
          enddo
        else
          do j=1,ny
            do i=1,nx
              Utau_t(i,j) = Utau_t_const
            enddo
          enddo
        endif
      elseif (Utau_t_flag.eq.1.0) then
        do j=1,ny
          do i=1,nx
            call surface_snow_1(tair_grid(i,j),windspd_grid(i,j),
     &        sprec(i,j),ro_soft_snow(i,j),Utau_t(i,j),
     &        ro_soft_snow_old(i,j),dt,snow_z0,ht_windobs,
     &        ro_nsnow(i,j))
          enddo
        enddo
        if (curve_lg_scale_flag.eq.1.0) then
          do j=1,ny
            do i=1,nx
              Utau_t(i,j) = curve_wt_lg(i,j) * Utau_t(i,j)
            enddo
          enddo
        endif
      endif

c Set the blowing snow flag to zero until it is clear that we will
c   have blowing snow.
      bs_flag = 0.0

c If the wind speed is lower that some threshold, then don't
c   need to to any of the snow transport computations.
      if (windspd_flag.ge.wind_min) then

c Get the wind direction indexing arrays for this particular
c   wind event (time step).
        call getdirection(nx,ny,uwind_grid,vwind_grid,index_ue,
     &    index_uw,index_vn,index_vs)

c Solve for Utau and z_0 if snow is saltating, else solve assuming
c   z_0 is known from snow depth and/or veg type, and solve for
c   Utau.
        call solveUtau(Utau,ht_windobs,windspd_grid,C_z,vonKarman,
     &    gravity,z_0,h_star,h_const,vegsnowd_xy,snow_d,
     &    snow_z0,veg_z0,bs_flag,nx,ny,Utau_t,soft_snow_d)

c If the blowing snow flag indicates wind transported snow
c   somewhere within the domain (bs_flag = 1.0), run the saltation
c   and suspension models.
        if (bs_flag.eq.1.0) then

c Solve for the saltation flux.
          print *,'         Saltation'
          call saltation(Qsalt,deltax,fetch,Utau,Utau_t,nx,ny,
     &      ro_air,gravity,vegsnowd_xy,snow_d,
     &      Qsalt_max,Qsalt_maxu,Qsalt_maxv,deltay,Qsalt_u,Qsalt_v,
     &      index_ue,index_uw,index_vn,index_vs,uwind_grid,
     &      vwind_grid,xmu,soft_snow_d,bc_flag)

c Solve for the suspension flux.
          print *,'         Suspension'
          call suspension(Utau,vonKarman,nx,ny,conc_salt,
     &      Qsalt,Qsusp,z_0,h_star,dz_susp,ztop_susp,pi,
     &      fall_vel,Ur_const,Up_const,Utau_t,Qsubl,ht_rhobs,
     &      tair_grid,rh_grid,Qsusp_u,Qsusp_v,uwind_grid,
     &      vwind_grid)

        elseif (bs_flag.eq.0.0) then

          call noblowsnow(nx,ny,Qsalt_max,Qsalt_maxu,
     &      Qsalt_maxv,Qsalt,Qsalt_u,Qsalt_v,dh_salt,dh_salt_u,
     &      dh_salt_v,conc_salt,Qsusp,Qsusp_u,Qsusp_v,dh_susp,
     &      dh_susp_u,dh_susp_v,Qsubl,dh_subgrid)

        endif

      else

c This 'noblowsnow' call zeros out data from a previous time step
c   that had blowing snow.
        call noblowsnow(nx,ny,Qsalt_max,Qsalt_maxu,
     &    Qsalt_maxv,Qsalt,Qsalt_u,Qsalt_v,dh_salt,dh_salt_u,
     &    dh_salt_v,conc_salt,Qsusp,Qsusp_u,Qsusp_v,dh_susp,
     &    dh_susp_u,dh_susp_v,Qsubl,dh_subgrid)

      endif

c Compute the new snow depth due to accumulation from precipitation,
c   saltation, and suspension, and the mass loss due to
c   sublimation.
      call accum(snow_d,nx,ny,ro_snow,dt,ro_water,
     &  deltax,deltay,vegtype,vegsnowd_xy,
     &  index_ue,index_uw,index_vn,index_vs,
     &  Qsalt_u,Qsalt_v,Qsusp_u,Qsusp_v,Qsubl,dh_salt,
     &  dh_salt_u,dh_salt_v,dh_susp,dh_susp_u,dh_susp_v,
     &  wbal_qsubl,wbal_salt,wbal_susp,bs_flag,
     &  soft_snow_d,topo,topo_land,topoflag,subgrid_flag,
     &  winddir_grid,tabler_nn,tabler_ss,tabler_ee,tabler_ww,
     &  tabler_ne,tabler_se,tabler_sw,tabler_nw,
     &  uwind_grid,vwind_grid,wbal_subgrid,sum_qsubl,
     &  sum_trans,swe_depth,snow_depth,ro_snow_grid,
     &  dh_subgrid,tabler_dir)

c Use the changes in swe due to saltation, suspension, and
c   blowing snow sublimation to adjust the multilayer snowpack
c   layers.
      if (multilayer_snowpack.eq.1) then
        do j=1,ny
          do i=1,nx
            swe_change = wbal_qsubl(i,j) + wbal_salt(i,j) +
     &        wbal_susp(i,j) + wbal_subgrid(i,j)

c Net mass loss for this grid cell at this time step.
            if (swe_change.lt.0.0) then
              swe_change_tmp = -swe_change

c Extract the vertical column for this i,j point, and send it
c   to the subroutine. *** Note that I should use f95, then I would

c   not have to do this (I could pass in subsections of the arrays).
              do k=1,nz_max
                swe_lyr_z(k) = swe_lyr(i,j,k)
              enddo

c Check to see whether a layer reduction is required.
              CALL REDUCE_LAYERS(swe_change_tmp,swe_lyr_z,JJ(i,j))

c Re-build the 3-D array.  See note above about using f95 to avoid this.
              do k=1,nz_max
                swe_lyr(i,j,k) = swe_lyr_z(k)
              enddo

c Update the snow layer thicknesses, and recalculate the total
c   snow and swe depths.  Assume this swe change does not change
c   the snow density and does not change the soft snow depth.  It
c   only reduces the snow depth and the associated swe depth.
              snow_depth(i,j) = 0.0
              swe_depth(i,j) = 0.0
              do k=1,JJ(i,j)
                dy_snow(i,j,k) = swe_lyr(i,j,k) * ro_water /
     &            ro_layer(i,j,k)
c               ro_layer(i,j,k) = ro_layer(i,j,k)
                snow_depth(i,j) = snow_depth(i,j) + dy_snow(i,j,k)
                swe_depth(i,j) = swe_depth(i,j) + swe_lyr(i,j,k)
              enddo

c Net mass gain for this grid cell at this time step.
            elseif (swe_change.gt.0.0) then

c Add to the existing top layer.
              swe_lyr(i,j,JJ(i,j)) = swe_lyr(i,j,JJ(i,j)) +
     &          swe_change
c             ro_layer(i,j,k) = ro_layer(i,j,k)
              dy_snow(i,j,JJ(i,j)) = swe_lyr(i,j,JJ(i,j)) *
     &          ro_water / ro_layer(i,j,JJ(i,j))

c Update the snow layer thicknesses, and recalculate the total
c   snow and swe depths.  Assume this swe change does not change
c   the snow density and does not change the soft snow depth.  It
c   only reduces the snow depth and the associated swe depth.
              snow_depth(i,j) = 0.0
              swe_depth(i,j) = 0.0
              do k=1,JJ(i,j)
                snow_depth(i,j) = snow_depth(i,j) + dy_snow(i,j,k)
                swe_depth(i,j) = swe_depth(i,j) + swe_lyr(i,j,k)
              enddo

            else

              snow_depth(i,j) = 0.0
              swe_depth(i,j) = 0.0
              do k=1,JJ(i,j)
                snow_depth(i,j) = snow_depth(i,j) + dy_snow(i,j,k)
                swe_depth(i,j) = swe_depth(i,j) + swe_lyr(i,j,k)
              enddo

            endif

          enddo
        enddo
      endif

c Perform a water balance check (see notes in this subroutine).
      if (seaice_run.eq.0.0) then
        call waterbal_snowtran(w_balance,prec,Qcs,
     &    runoff,d_canopy_int,swe_depth,glacier_melt,iter,
     &    wbal_qsubl,wbal_salt,wbal_susp,wbal_subgrid,nx,ny,
     &    swe_depth_old,swesublim,canopy_unload,canopy_int,
     &    canopy_int_old)

c       call waterbal_snowtran_sums(w_balance,sum_prec,sum_Qcs,
c    &    sum_runoff,canopy_int,swe_depth,sum_glacmelt,iter,
c    &    sum_qsubl,sum_trans,nx,ny,snow_d_init,ro_snow,ro_water,
c    &    sum_sfcsublim)
      endif

c If this is a sea ice run, zero out the ocean grid cells that
c   have no sea ice in them.
      if (seaice_run.ne.0.0) then
        CALL ZERO_SEAICE_SNOW(nx,ny,snow_depth,ro_snow_grid,
     &    ro_snow,swe_depth,swe_depth_old,canopy_int_old,JJ,
     &    tslsnowfall,dy_snow,swe_lyr,ro_layer,T_old,
     &    multilayer_snowpack,tsls_threshold,seaice_conc)
      endif

c Save the mass balance variables from this time step.
      do j=1,ny
        do i=1,nx
          swe_depth_old(i,j) = swe_depth(i,j)
          canopy_int_old(i,j) = canopy_int(i,j)
        enddo
      enddo

c In SnowTran-3D, the summed snow precipitation were in units
c   of snow-depth.  The rest of the routines assume that it is in
c   swe units.
      do j=1,ny
        do i=1,nx
          sum_sprec(i,j) = sum_sprec(i,j) * ro_snow / ro_water
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine waterbal_snowtran(w_balance,prec,Qcs,
     &  runoff,d_canopy_int,swe_depth,glacier_melt,iter,
     &  wbal_qsubl,wbal_salt,wbal_susp,wbal_subgrid,nx,ny,
     &  swe_depth_old,swesublim,canopy_unload,canopy_int,
     &  canopy_int_old)

      implicit none

      include 'snowmodel.inc'

      integer iter,nx,ny,i,j

      real w_balance(nx_max,ny_max),prec(nx_max,ny_max),
     &  Qcs(nx_max,ny_max),runoff(nx_max,ny_max),
     &  d_canopy_int(nx_max,ny_max),swe_depth(nx_max,ny_max),
     &  glacier_melt(nx_max,ny_max),wbal_qsubl(nx_max,ny_max),
     &  wbal_salt(nx_max,ny_max),swe_depth_old(nx_max,ny_max),
     &  swesublim(nx_max,ny_max),wbal_susp(nx_max,ny_max),
     &  wbal_subgrid(nx_max,ny_max),canopy_unload(nx_max,ny_max),
     &  canopy_int_old(nx_max,ny_max),canopy_int(nx_max,ny_max)

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

c The subroutine WATERBAL_SNOWTRAN is used if the model simulation
c   includes SnowTran-3D.
      do j=1,ny
        do i=1,nx
          w_balance(i,j) = swe_depth_old(i,j) - swe_depth(i,j) +
     &      prec(i,j) - runoff(i,j) + glacier_melt(i,j) +
     &      wbal_qsubl(i,j) + wbal_salt(i,j) + wbal_susp(i,j) +
     &      wbal_subgrid(i,j) - swesublim(i,j) + canopy_int_old(i,j) -
     &      canopy_int(i,j) + Qcs(i,j)

          if (abs(w_balance(i,j)).gt.1.0e-5)
     &      print*,'water imbalance at iter =',iter,' ',w_balance(i,j)

        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine waterbal_snowtran_sums(w_balance,sum_prec,sum_Qcs,
     &  sum_runoff,canopy_int,swe_depth,sum_glacmelt,iter,
     &  sum_qsubl,sum_trans,nx,ny,snow_d_init,ro_snow,ro_water,
     &  sum_sfcsublim)

      implicit none

      include 'snowmodel.inc'

      integer iter,nx,ny,i,j

      real w_balance(nx_max,ny_max),sum_prec(nx_max,ny_max),
     &  sum_Qcs(nx_max,ny_max),sum_runoff(nx_max,ny_max),
     &  canopy_int(nx_max,ny_max),swe_depth(nx_max,ny_max),
     &  sum_glacmelt(nx_max,ny_max),sum_qsubl(nx_max,ny_max),
     &  sum_trans(nx_max,ny_max),snow_d_init(nx_max,ny_max),
     &  sum_sfcsublim(nx_max,ny_max)

      real ro_snow,ro_water

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

c The subroutine WATERBAL_SNOWTRAN is used if the model simulation
c   includes SnowTran-3D.
      do j=1,ny
        do i=1,nx
          w_balance(i,j) = sum_prec(i,j) + sum_Qcs(i,j) - 
     &      sum_runoff(i,j) - canopy_int(i,j) - swe_depth(i,j) +
     &      sum_glacmelt(i,j) + sum_qsubl(i,j) + sum_trans(i,j) +
     &      snow_d_init(i,j) * ro_snow/ro_water -
     &      sum_sfcsublim(i,j)

          if (abs(w_balance(i,j)).gt.1.0e-4)
     &      print*,'water imbalance at iter =',iter,' ',w_balance(i,j)
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine noblowsnow(nx,ny,Qsalt_max,Qsalt_maxu,
     &  Qsalt_maxv,Qsalt,Qsalt_u,Qsalt_v,dh_salt,dh_salt_u,
     &  dh_salt_v,conc_salt,Qsusp,Qsusp_u,Qsusp_v,dh_susp,
     &  dh_susp_u,dh_susp_v,Qsubl,dh_subgrid)

      implicit none

      include 'snowmodel.inc'

      integer nx,ny,i,j
      real Qsalt_max(nx_max,ny_max)
      real Qsalt_maxu(nx_max,ny_max),Qsalt_maxv(nx_max,ny_max)
      real Qsalt(nx_max,ny_max)
      real Qsalt_u(nx_max,ny_max),Qsalt_v(nx_max,ny_max)
      real dh_salt(nx_max,ny_max)
      real dh_salt_u(nx_max,ny_max),dh_salt_v(nx_max,ny_max)

      real conc_salt(nx_max,ny_max)

      real Qsusp(nx_max,ny_max)
      real Qsusp_u(nx_max,ny_max),Qsusp_v(nx_max,ny_max)
      real dh_susp(nx_max,ny_max)
      real dh_susp_u(nx_max,ny_max),dh_susp_v(nx_max,ny_max)
      real dh_subgrid(nx_max,ny_max)

      real Qsubl(nx_max,ny_max)

      do i=1,nx
        do j=1,ny
          Qsalt_max(i,j) = 0.0
          Qsalt_maxu(i,j) = 0.0
          Qsalt_maxv(i,j) = 0.0
          Qsalt(i,j) = 0.0
          Qsalt_u(i,j) = 0.0
          Qsalt_v(i,j) = 0.0
          dh_salt(i,j) = 0.0
          dh_salt_u(i,j) = 0.0
          dh_salt_v(i,j) = 0.0
          conc_salt(i,j) = 0.0
          Qsusp(i,j) = 0.0
          Qsusp_u(i,j) = 0.0
          Qsusp_v(i,j) = 0.0
          dh_susp(i,j) = 0.0
          dh_susp_u(i,j) = 0.0
          dh_susp_v(i,j) = 0.0
          Qsubl(i,j) = 0.0
          dh_subgrid(i,j) = 0.0
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine precip(snow_d,sprec,nx,ny,ro_snow,
     &  ro_water,sum_sprec,soft_snow_d,sum_prec)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny

      real ro_snow,ro_water

      real sprec(nx_max,ny_max)
      real snow_d(nx_max,ny_max)
      real sum_sprec(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)
      real sum_prec(nx_max,ny_max)

      do i=1,nx
        do j=1,ny

c Place the new snow in the soft snow layer.
          soft_snow_d(i,j) = soft_snow_d(i,j) +
     &      sprec(i,j) * ro_water / ro_snow

c Update the snow depth resulting from swe precipitation.
          snow_d(i,j) = snow_d(i,j) + sprec(i,j) * ro_water / ro_snow

c Sum the precipitation in terms of snow depth.
          sum_sprec(i,j) = sum_sprec(i,j) + sprec(i,j) *
     &      ro_water / ro_snow

c Sum the precipitation in terms of swe.
          sum_prec(i,j) = sum_prec(i,j) + sprec(i,j)

        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine twolayer1(nx,ny,soft_snow_d,tair_grid)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny

      real tmax,tfreeze

      real tair_grid(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)

c Place the soft snow layer in the hard snow layer if the air
c   temperature gets too high.
      tmax = 3.0
      tfreeze = 273.16
        do i=1,nx
          do j=1,ny
            if (tair_grid(i,j).ge.tfreeze+tmax) soft_snow_d(i,j) = 0.0
          enddo
        enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine accum(snow_d,nx,ny,ro_snow,dt,ro_water,
     &  deltax,deltay,vegtype,vegsnowd_xy,
     &  index_ue,index_uw,index_vn,index_vs,
     &  Qsalt_u,Qsalt_v,Qsusp_u,Qsusp_v,Qsubl,dh_salt,
     &  dh_salt_u,dh_salt_v,dh_susp,dh_susp_u,dh_susp_v,
     &  wbal_qsubl,wbal_salt,wbal_susp,bs_flag,
     &  soft_snow_d,topo,topo_land,topoflag,subgrid_flag,
     &  winddir_grid,tabler_nn,tabler_ss,tabler_ee,tabler_ww,
     &  tabler_ne,tabler_se,tabler_sw,tabler_nw,
     &  uwind_grid,vwind_grid,wbal_subgrid,sum_qsubl,
     &  sum_trans,swe_depth,snow_depth,ro_snow_grid,
     &  dh_subgrid,tabler_dir)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny

      real ro_snow,dt,deltax,deltay,bs_flag,topoflag,ro_water
      real snowdmin,hard_snow_d,subgrid_flag,tabler_dir

      real snow_d(nx_max,ny_max)
      real snow_d_tmp(nx_max,ny_max)
      real snow_depth(nx_max,ny_max)
      real swe_depth(nx_max,ny_max)
      real ro_snow_grid(nx_max,ny_max)
      real snow_d_tabler(nx_max,ny_max)

      real tabler_nn(nx_max,ny_max)
      real tabler_ss(nx_max,ny_max)
      real tabler_ee(nx_max,ny_max)
      real tabler_ww(nx_max,ny_max)
      real tabler_ne(nx_max,ny_max)
      real tabler_se(nx_max,ny_max)
      real tabler_sw(nx_max,ny_max)
      real tabler_nw(nx_max,ny_max)
      real winddir_grid(nx_max,ny_max)
      real uwind_grid(nx_max,ny_max)
      real vwind_grid(nx_max,ny_max)

      real soft_snow_d(nx_max,ny_max)
      real Qsubl(nx_max,ny_max)
      real topo(nx_max,ny_max)
      real topo_land(nx_max,ny_max)

      real dh_salt(nx_max,ny_max)
      real dh_salt_u(nx_max,ny_max)
      real dh_salt_v(nx_max,ny_max)

      real dh_susp(nx_max,ny_max)
      real dh_susp_u(nx_max,ny_max)
      real dh_susp_v(nx_max,ny_max)

      real dh_subgrid(nx_max,ny_max)

      real Qsalt_u(nx_max,ny_max)
      real Qsalt_v(nx_max,ny_max)

      real Qsusp_u(nx_max,ny_max)
      real Qsusp_v(nx_max,ny_max)

      real wbal_qsubl(nx_max,ny_max)
      real wbal_salt(nx_max,ny_max)
      real wbal_susp(nx_max,ny_max)
      real wbal_subgrid(nx_max,ny_max)
      real sum_qsubl(nx_max,ny_max)
      real sum_trans(nx_max,ny_max)

      real vegtype(nx_max,ny_max)
      real vegsnowd_xy(nx_max,ny_max)

      integer index_ue(ny_max,2*nx_max+1)
      integer index_uw(ny_max,2*nx_max+1)
      integer index_vn(nx_max,2*ny_max+1)
      integer index_vs(nx_max,2*ny_max+1)

c COMPUTE THE NEW SNOW DEPTH.

c PRECIPITATION
c Account for the addition due to snow precipitation.
c This is now updated at the beginning of the program (day).

c Sum the precipitation in terms of snow depth.
      if (bs_flag.eq.1.0) then

c SALTATION

      call getnewdepth(nx,ny,deltax,deltay,Qsalt_u,
     &  Qsalt_v,dh_salt_u,dh_salt_v,index_ue,index_uw,
     &  index_vn,index_vs,ro_snow,dt,vegsnowd_xy,snow_d,
     &  soft_snow_d)

      do i=1,nx
        do j=1,ny
          dh_salt(i,j) = dh_salt_u(i,j) + dh_salt_v(i,j)
        enddo
      enddo

c SUSPENSION

      call getnewdepth(nx,ny,deltax,deltay,Qsusp_u,
     &  Qsusp_v,dh_susp_u,dh_susp_v,index_ue,index_uw,
     &  index_vn,index_vs,ro_snow,dt,vegsnowd_xy,snow_d,
     &  soft_snow_d)

      do i=1,nx
        do j=1,ny
          dh_susp(i,j) = dh_susp_u(i,j) + dh_susp_v(i,j)
        enddo
      enddo

c SUBLIMATION
c Make adjustments for the case where there is no snow available
c   on the ground (or captured within the vegetation) to be
c   eroded.
      do i=1,nx
        do j=1,ny
          hard_snow_d = snow_d(i,j) - soft_snow_d(i,j)
          snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)

c Convert Qsubl to sublimated snow depth.
          Qsubl(i,j) = Qsubl(i,j) * dt / ro_snow

          if (snow_d(i,j).gt.snowdmin) then
            if (snow_d(i,j)+Qsubl(i,j).le.snowdmin) then
              Qsubl(i,j) = snowdmin - snow_d(i,j)
            endif
          else
            Qsubl(i,j) = 0.0
          endif
        enddo
      enddo

c Save a copy of the snow distribution to be used to calculate the
c   snow distribution changes resulting from the subgrid
c   redistribution.
      do i=1,nx
        do j=1,ny
          snow_d_tmp(i,j) = snow_d(i,j)
        enddo
      enddo

c Run the subgrid parameterization to account for unrealistic
c   snow accumulation spikes.
      if (subgrid_flag.eq.1.0) then

c Do the Tabler corrections while considering the SnowPack density
c   contribution to snow depth.  When done, convert back to SnowTran
c   constant density convention.
        do i=1,nx
          do j=1,ny
            snow_d_tabler(i,j) = snow_d(i,j) * ro_snow /
     &        ro_snow_grid(i,j)
          enddo
        enddo

        call subgrid_1(nx,ny,snow_d_tabler,
     &    index_ue,index_uw,index_vn,index_vs,
     &    tabler_nn,tabler_ss,tabler_ee,tabler_ww,
     &    tabler_ne,tabler_se,tabler_sw,tabler_nw,uwind_grid,
     &    vwind_grid,tabler_dir)

        do i=1,nx
          do j=1,ny
            snow_d(i,j) = snow_d_tabler(i,j) * ro_snow_grid(i,j) /
     &        ro_snow
          enddo
        enddo

      endif

c Calculate the snow depth resulting from the subgrid
c   redistribution.
      do i=1,nx
        do j=1,ny
          dh_subgrid(i,j) = snow_d(i,j) - snow_d_tmp(i,j)
        enddo
      enddo

c Account for decreases in snow depth due to sublimation.
      do i=1,nx
        do j=1,ny
          snow_d(i,j) = snow_d(i,j) + Qsubl(i,j)
          soft_snow_d(i,j) = soft_snow_d(i,j) + Qsubl(i,j)
        enddo
      enddo

c Update the surface topography resulting from the snow setting
c   on the land.
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

c MOISTURE BALANCE
c Save enough information to do a moisture balance.
      do i=1,nx
        do j=1,ny

c Save the sublimation in terms of swe depth.
          wbal_qsubl(i,j) = Qsubl(i,j) * ro_snow / ro_water

c Save the saltation in terms of swe depth.
          wbal_salt(i,j) = dh_salt(i,j) * ro_snow / ro_water

c Save the suspension in terms of swe depth.
          wbal_susp(i,j) = dh_susp(i,j) * ro_snow / ro_water

c Save the subgrid redistribution in terms of swe depth.
          wbal_subgrid(i,j) = dh_subgrid(i,j) * ro_snow / ro_water

c Fill summing arrays of the sublimation and transport quantities.
          sum_qsubl(i,j) = sum_qsubl(i,j) + wbal_qsubl(i,j)
          sum_trans(i,j) = sum_trans(i,j) + wbal_salt(i,j) +
     &      wbal_susp(i,j) + wbal_subgrid(i,j)

        enddo
      enddo

      do i=1,nx
        do j=1,ny

c Convert any snow-depth adjustments that occurred in SnowTran-3D
c   to swe (using the SnowTran-3D constant snow density) so
c   that it can be used in SNOWPACK that accounts for the
c   time-evolution of snow density.
          swe_depth(i,j) = snow_d(i,j) * ro_snow / ro_water

c Calculate the snow depth using the spatially-distributed snow density
c   from the snowpack model.
          snow_depth(i,j) = swe_depth(i,j) *
     &      ro_water / ro_snow_grid(i,j)

        enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine smoother9(nx,ny,snow)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny
      real snow(nx_max,ny_max)
      real snow_tmp(nx_max,ny_max)

c Performs a 9-point smoothing operation.

c The result at each grid point is a weighted average of the grid
c   point and the surrounding 8 points.  The center point receives
c   a weight of 1.0, the points at each side and above and below
c   receive a weight of 0.5, and corner points receive a weight of
c   0.3.  All points are multiplied by their weights and summed,
c   then divided by the total weight.

c Do the interior.
      do i=2,nx-1
        do j=2,ny-1
          snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) +
     &      snow(i,j+1) + snow(i-1,j) + snow(i+1,j)) + 0.3 *
     &      (snow(i-1,j-1) + snow(i+1,j+1) + snow(i-1,j+1) +
     &      snow(i+1,j-1))) / 4.2
        enddo
      enddo

c Do the sides.
      j = 1
      do i=2,nx-1
        snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i-1,j) +
     &    snow(i+1,j)) + 0.3 * (snow(i+1,j+1) + snow(i-1,j+1))) / 3.1
      enddo

      j = ny
      do i=2,nx-1
        snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i-1,j) +
     &    snow(i+1,j)) + 0.3 * (snow(i+1,j-1) + snow(i-1,j-1))) / 3.1
      enddo

      i = 1
      do j=2,ny-1
        snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i,j+1) +
     &    snow(i+1,j)) + 0.3 * (snow(i+1,j-1) + snow(i+1,j+1))) / 3.1
      enddo

      i = nx
      do j=2,ny-1
        snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i,j+1) +
     &    snow(i-1,j)) + 0.3 * (snow(i-1,j-1) + snow(i-1,j+1))) / 3.1
      enddo

c Do the corners.
      i = 1
      j = 1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i+1,j)) +
     &  0.3 * snow(i+1,j+1)) / 2.3

      i = nx
      j = 1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i-1,j)) +
     &  0.3 * snow(i-1,j+1)) / 2.3

      i = 1
      j = ny
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i+1,j)) +
     &  0.3 * snow(i+1,j-1)) / 2.3

      i = nx
      j = ny
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i-1,j)) +
     &  0.3 * snow(i-1,j-1)) / 2.3

c Return the smoothed array.
      do i=1,nx
        do j=1,ny
          snow(i,j) = snow_tmp(i,j)
        enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine suspension(Utau,vonKarman,nx,ny,conc_salt,
     &  Qsalt,Qsusp,z_0,h_star,dz,ztop,pi,
     &  fall_vel,Ur_const,Up_const,Utau_t,Qsubl,ht_rhobs,
     &  tair_grid,rh_grid,Qsusp_u,Qsusp_v,uwind_grid,
     &  vwind_grid)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,nzsteps,iz

      real vonKarman,dz,ztop,fall_vel,Ur_const,Up_const
      real ht_rhobs,V_susp,V_salt,pi
      real U_p,Utau_fallvel,U_r,phistar_Cr,product,conc,z

      real Utau(nx_max,ny_max)
      real Utau_t(nx_max,ny_max)
      real uwind_grid(nx_max,ny_max)
      real vwind_grid(nx_max,ny_max)
      real z_0(nx_max,ny_max)
      real h_star(nx_max,ny_max)
      real conc_salt(nx_max,ny_max)
      real Qsalt(nx_max,ny_max)
      real Qsusp(nx_max,ny_max)
      real Qsusp_u(nx_max,ny_max)
      real Qsusp_v(nx_max,ny_max)
      real Qsubl(nx_max,ny_max)

      real tair_grid(nx_max,ny_max)
      real rh_grid(nx_max,ny_max)

c Compute the mass concentration of suspended snow according to
c   Kind (1992).

      do i=1,nx
      do j=1,ny
        if (Qsalt(i,j).gt.0.0) then
          Utau_fallvel = Utau(i,j) / fall_vel
          if (h_star(i,j).eq.z_0(i,j)) h_star(i,j) = 2.0 * z_0(i,j)
          U_r = Utau(i,j)/vonKarman * log(h_star(i,j)/z_0(i,j))
          phistar_Cr = Utau(i,j)/U_r * Ur_const
          product = phistar_Cr * Utau_fallvel
          U_p = Up_const * Utau_t(i,j)

c Compute the concentration in the saltation layer (kg/m**3).
          conc_salt(i,j) = Qsalt(i,j) / (h_star(i,j) * U_p)

          nzsteps = int((ztop - h_star(i,j)) / dz)

          Qsusp(i,j) = 0.0
          Qsubl(i,j) = 0.0

          do iz=1,nzsteps
            z = h_star(i,j) + 0.5 * dz + real(iz - 1) * dz

c Compute the concentration of the suspended snow at height z.
            conc = conc_salt(i,j) * ((product + 1.0) *
     &        (z/h_star(i,j))**((-fall_vel)/(vonKarman*Utau(i,j))) -
     &        product)
            conc = max(conc,0.0)

c Only do The integration if the concentration is non-zero.
            if (conc.gt.0.0) then

c Compute the sublimation due to suspension.
              call getsublim(z,rh_grid(i,j),tair_grid(i,j),Utau(i,j),
     &          z_0(i,j),V_susp,V_salt,Utau_t(i,j),ht_rhobs,1.0,pi)

c Perform the quadrature (summation), without the constants.
              if (z.eq.z_0(i,j)) z = 1.2 * z_0(i,j)
              Qsusp(i,j) = Qsusp(i,j) + conc * log(z/z_0(i,j)) * dz
              Qsubl(i,j) = Qsubl(i,j) + conc * V_susp * dz

            endif

          enddo

c Finish the quadratures.
c Include the constants for Qsusp.
        Qsusp(i,j) = Utau(i,j) / vonKarman * Qsusp(i,j)

c Include the sublimation contribution due to saltation.
        z = h_star(i,j) / 2.0
        call getsublim(z,rh_grid(i,j),tair_grid(i,j),Utau(i,j),
     &    z_0(i,j),V_susp,V_salt,Utau_t(i,j),ht_rhobs,0.0,pi)

        Qsubl(i,j) = Qsubl(i,j) +
     &    V_salt * conc_salt(i,j) * h_star(i,j)

        else
          conc_salt(i,j) = 0.0
          Qsusp(i,j) = 0.0
          Qsubl(i,j) = 0.0
        endif

      enddo
      enddo

c Separate the east-west and the north-south suspended transport
c   components; the vector sum should equal Qsusp.
      do i=1,nx
        do j=1,ny
          Qsusp_u(i,j) = Qsusp(i,j) * abs(uwind_grid(i,j)) /
     &      sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
          Qsusp_v(i,j) = Qsusp(i,j) * abs(vwind_grid(i,j)) /
     &      sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
        enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine saltation(Qsalt,deltax,fetch,Utau,Utau_t,nx,ny,
     &  ro_air,gravity,vegsnowd_xy,snow_d,
     &  Qsalt_max,Qsalt_maxu,Qsalt_maxv,deltay,Qsalt_u,Qsalt_v,
     &  index_ue,index_uw,index_vn,index_vs,uwind_grid,
     &  vwind_grid,xmu,soft_snow_d,bc_flag)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny
      integer k,istart,iend,jstart,jend

      real deltax,deltay,fetch,ro_air,gravity,dUtau,xmu,
     &  blowby,bc_flag

      real Qsalt_max(nx_max,ny_max)
      real Qsalt_maxu(nx_max,ny_max)
      real Qsalt_maxv(nx_max,ny_max)
      real Qsalt(nx_max,ny_max)
      real Qsalt_u(nx_max,ny_max)
      real Qsalt_v(nx_max,ny_max)
      real Utau(nx_max,ny_max)
      real Utau_t(nx_max,ny_max)
      real uwind_grid(nx_max,ny_max)
      real vwind_grid(nx_max,ny_max)
      real snow_d(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)
      real vegsnowd_xy(nx_max,ny_max)

      integer index_ue(ny_max,2*nx_max+1)
      integer index_uw(ny_max,2*nx_max+1)
      integer index_vn(nx_max,2*ny_max+1)
      integer index_vs(nx_max,2*ny_max+1)

      real scale_EW,scale_NS

c The blowby parameter is implemented to account for the erosion
c   of the tops of deep snow accumulations.  It corrects a
c   deficiency in the du*/dx* < 0 formulation.  It is a number that
c   should range from 0 to 1.0, and represents the fraction of the
c   upwind saltation flux that is transfered farther downwind into
c   the next grid cell.  So, the bigger the number, the less
c   peaked the drift accumulation profile is.  blowby = 0.0 is the
c   original model.  I am now using the Tabler surfaces to do the
c   same kind of thing, so here I hard-code the parameter as in the
c   original model.
      blowby = 0.0

c Compute the maximum possible saltation flux, assuming that
c   an abundance of snow is available at the surface.
      do i=1,nx
        do j=1,ny

c For a given wind speed, find Qsalt_max.
          Qsalt_max(i,j) = 0.68 * ro_air / gravity *
     &      Utau_t(i,j) / Utau(i,j) * (Utau(i,j)**2 - Utau_t(i,j)**2)
          Qsalt_max(i,j) = max(Qsalt_max(i,j),0.0)

c Now weight the max saltation flux for the u and v wind
c   components, where the vector sum should equal Qsalt_max.
          Qsalt_maxu(i,j) = Qsalt_max(i,j) * abs(uwind_grid(i,j)) /
     &      sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
          Qsalt_maxv(i,j) = Qsalt_max(i,j) * abs(vwind_grid(i,j)) /
     &      sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
        enddo
      enddo

c Define an upwind boundary condition.  If bc_flag = 1.0 then it is
c   assumed that the inflow saltation flux has reached steady state.
c   If bc_flag = 0.0 then the saltation flux is assumed to be zero.
c   The boundary condition is implemented by initializing the arrays
c   to Qsalt_max, and since upwind boundaries are not called in
c   the Qsalt computation, they stay in effect for the future 
c   accumulation/erosion computation.
      if (bc_flag.eq.0.0) then
        do i=1,nx
          do j=1,ny
c Zero incoming flux at the boundaries.
            Qsalt_u(i,j) = 0.0
            Qsalt_v(i,j) = 0.0
          enddo
        enddo
      elseif (bc_flag.eq.1.0) then
        do i=1,nx
          do j=1,ny
c Steady-state (maximum) incoming flux at the boundaries.
            Qsalt_u(i,j) = Qsalt_maxu(i,j)
            Qsalt_v(i,j) = Qsalt_maxv(i,j)
          enddo
        enddo
      endif

c Define the scaling coefficients for Eqn. 9 in L&S 1998. Don't
c   let them be greater than 1.0 or you will make more snow than
c   there was before.
      scale_EW =  xmu * deltax / fetch
      scale_EW = min(1.0,scale_EW)
      scale_NS =  xmu * deltay / fetch
      scale_NS = min(1.0,scale_NS)

c Consider WESTERLY winds.
      do j=1,ny
        do k=1,index_uw(j,1)
          istart = index_uw(j,k*2)+1
          iend = index_uw(j,k*2+1)
          do i=istart,iend
            dUtau = Utau(i,j) - Utau(i-1,j)
            if (dUtau.ge.0.0) then
              Qsalt_u(i,j) = Qsalt_u(i-1,j) + scale_EW *
     &          (Qsalt_maxu(i,j) - Qsalt_u(i-1,j))
            else
c             Qsalt_u(i,j) = min(Qsalt_u(i-1,j),Qsalt_maxu(i,j))

              if (Qsalt_u(i-1,j).lt.Qsalt_maxu(i,j)) then
                Qsalt_u(i,j) = Qsalt_u(i-1,j)
              else
                Qsalt_u(i,j) =
     &            max(blowby*Qsalt_u(i-1,j),Qsalt_maxu(i,j))
              endif

            endif
          enddo
        enddo
      enddo

c Consider EASTERLY winds.
      do j=1,ny
        do k=1,index_ue(j,1)
          iend = index_ue(j,k*2)
          istart = index_ue(j,k*2+1)-1
          do i=istart,iend,-1
            dUtau = Utau(i,j) - Utau(i+1,j)
            if (dUtau.ge.0.0) then
              Qsalt_u(i,j) = Qsalt_u(i+1,j) + scale_EW *
     &          (Qsalt_maxu(i,j) - Qsalt_u(i+1,j))
            else
c             Qsalt_u(i,j) = min(Qsalt_u(i+1,j),Qsalt_maxu(i,j))

              if (Qsalt_u(i+1,j).lt.Qsalt_maxu(i,j)) then
                Qsalt_u(i,j) = Qsalt_u(i+1,j)
              else
                Qsalt_u(i,j) =
     &            max(blowby*Qsalt_u(i+1,j),Qsalt_maxu(i,j))
              endif

            endif
          enddo
        enddo
      enddo

c Consider SOUTHERLY winds.
      do i=1,nx
        do k=1,index_vs(i,1)
          jstart = index_vs(i,k*2)+1
          jend = index_vs(i,k*2+1)
          do j=jstart,jend
            dUtau = Utau(i,j) - Utau(i,j-1)
            if (dUtau.ge.0.0) then
              Qsalt_v(i,j) = Qsalt_v(i,j-1) + scale_NS *
     &          (Qsalt_maxv(i,j) - Qsalt_v(i,j-1))
            else
c             Qsalt_v(i,j) = min(Qsalt_v(i,j-1),Qsalt_maxv(i,j))

              if (Qsalt_v(i,j-1).lt.Qsalt_maxv(i,j)) then
                Qsalt_v(i,j) = Qsalt_v(i,j-1)
              else
                Qsalt_v(i,j) =
     &            max(blowby*Qsalt_v(i,j-1),Qsalt_maxv(i,j))
              endif

            endif
          enddo
        enddo
      enddo

c Consider NORTHERLY winds.
      do i=1,nx
        do k=1,index_vn(i,1)
          jend = index_vn(i,k*2)
          jstart = index_vn(i,k*2+1)-1
          do j=jstart,jend,-1
            dUtau = Utau(i,j) - Utau(i,j+1)
            if (dUtau.ge.0.0) then
              Qsalt_v(i,j) = Qsalt_v(i,j+1) + scale_NS *
     &          (Qsalt_maxv(i,j) - Qsalt_v(i,j+1))
            else
c             Qsalt_v(i,j) = min(Qsalt_v(i,j+1),Qsalt_maxv(i,j))

              if (Qsalt_v(i,j+1).lt.Qsalt_maxv(i,j)) then
                Qsalt_v(i,j) = Qsalt_v(i,j+1)
              else
                Qsalt_v(i,j) =
     &            max(blowby*Qsalt_v(i,j+1),Qsalt_maxv(i,j))
              endif

            endif
          enddo
        enddo
      enddo

c Combine the u and v components to yield the total saltation flux
c   at each grid cell.
      do i=1,nx
        do j=1,ny
          Qsalt(i,j) = Qsalt_u(i,j) + Qsalt_v(i,j)
        enddo
      enddo

c Adjust Qsalt to account for the availablity of snow for transport;
c   taking into consideration whether there is snow on the ground,
c   the holding depth of the vegetation, etc..
      do i=1,nx
        do j=1,ny
          if (snow_d(i,j).le.vegsnowd_xy(i,j)) Qsalt(i,j) = 0.0
          if (soft_snow_d(i,j).le.0.0) Qsalt(i,j) = 0.0
        enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine solveUtau(Utau,ht_windobs,windspd_grid,C_z,vonKarman,
     &  gravity,z_0,h_star,h_const,vegsnowd_xy,snow_d,
     &  snow_z0,veg_z0,bs_flag,nx,ny,Utau_t,soft_snow_d)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny

      real bs_flag,guess,sfrac,vonKarman,ht_windobs,C_z,gravity
      real h_const,snow_z0,Utautmp,windtmp,wind_max
      real threshold,threshold_flag,z_0_tmp

      real Utau(nx_max,ny_max)
      real Utau_t(nx_max,ny_max)
      real windspd_grid(nx_max,ny_max)
      real z_0(nx_max,ny_max)
      real h_star(nx_max,ny_max)
      real snow_d(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)
      real veg_z0(nx_max,ny_max)
      real vegsnowd_xy(nx_max,ny_max)

c Initially set the blowing snow flag to no blowing snow
c   (bs_flag = 0.0).  Then, if snow is found to blow in any
c   domain grid cell, set the flag to on (bs_flag = 1.0).
      bs_flag = 0.0

c Build the Utau array.
      guess = 0.1
      do i=1,nx
      do j=1,ny

c Determine whether snow is saltating (this influences how Utau
c   and z_0 are computed).
        if (snow_d(i,j).le.vegsnowd_xy(i,j)) then

c Saltation will not occur.
          sfrac = snow_d(i,j) / max(vegsnowd_xy(i,j),veg_z0(i,j))
          z_0(i,j) = sfrac * snow_z0 + (1.0 - sfrac) * veg_z0(i,j)
          z_0_tmp = min(0.25*ht_windobs,z_0(i,j))
          Utau(i,j) = windspd_grid(i,j) *
     &      vonKarman / log(ht_windobs/z_0_tmp)
          h_star(i,j) = z_0(i,j) * h_const / C_z
        elseif (soft_snow_d(i,j).le.0.0) then
c Saltation will not occur.
          z_0(i,j) = snow_z0
          Utau(i,j) = windspd_grid(i,j) *
     &      vonKarman / log(ht_windobs/z_0(i,j))
          h_star(i,j) = z_0(i,j)
        else
c Saltation may occur.  Test for that possibility by assuming that
c   saltation is present, solving for Utau and z_0, and comparing
c   whether Utau exceeds Utau_t.  If it does not, set z_0 to that
c   of snow and recompute Utau.

c To help insure that the iteration converges, set the minimum
c   wind speed to be 1.0 m/s, and the maximum wind speed to be
c   30 m/s at 10-m height.
          windtmp = max(1.0,windspd_grid(i,j))
          wind_max = 30.0 * log(ht_windobs/snow_z0)/log(10.0/snow_z0)
          windtmp = min(windtmp,wind_max) 

c For u* over 0.6, use the relation z0 = 0.00734 u* - 0.0022,
c   instead of Equation (5) in Liston and Sturm (1998).  Note that
c   for windspeeds greater than about 35 m/s this will have to be
c   modified for the solution algorithm to converge (because the
c   roughness length will start to be higher than the obs height!).
          threshold = 0.6/vonKarman * log(ht_windobs/0.0022)
          if (windtmp.le.threshold) then
            threshold_flag = 1.0
          else
            threshold_flag = 2.0
          endif

          call solve1(Utautmp,guess,ht_windobs,windtmp,C_z,vonKarman,
     &      gravity,threshold_flag)

          if (Utautmp.gt.Utau_t(i,j)) then

c We have saltation.
            Utau(i,j) = Utautmp
            z_0(i,j) = C_z * Utau(i,j)**2 / (2.0 * gravity)
            h_star(i,j) = h_const * Utau(i,j)**2 / (2.0 * gravity)
            bs_flag = 1.0
          else

c We do not have saltation, but the vegetation is covered by snow.
c   Because we have determined that we do not have saltation, make
c   sure Utau does not exceed Utau_t.
            z_0(i,j) = snow_z0
            Utau(i,j) = windspd_grid(i,j) *
     &        vonKarman / log(ht_windobs/z_0(i,j))
            Utau(i,j) = min(Utau(i,j),Utau_t(i,j))
            h_star(i,j) = z_0(i,j) * h_const / C_z

          endif
        endif

      enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine solve1(xnew,guess,z,windtmp,C_z,vonKarman,
     &  gravity,threshold_flag)

      implicit none

      integer i,maxiter

      real xnew,guess,z,windtmp,C_z,vonKarman,tol,old,gravity
      real fprime,funct,threshold_flag

      tol = 1.0e-3
      maxiter = 20
      old = guess

      if (threshold_flag.eq.1.0) then

        do i=1,maxiter
          fprime = - 1.0 + 2.0 / old * windtmp * vonKarman *
     &      (log(z) - log(C_z/(2.0*gravity)) - 2.0*log(old))**(-2)
          funct = - old + windtmp * vonKarman *
     &      (log(z) - log(C_z/(2.0*gravity)) - 2.0*log(old))**(-1)
          xnew = old - funct/fprime
          if (abs(xnew - old).lt.tol) return
          old = xnew
        end do

      elseif (threshold_flag.eq.2.0) then

        old = 0.6
        do i=1,maxiter
          fprime = - 1.0 + windtmp * vonKarman *
     &      0.00734 / (0.00734 * old - 0.0022) *
     &      (log(z) - log(0.00734 * old - 0.0022))**(-2)
          funct = - old + windtmp * vonKarman *
     &      (log(z) - log(0.00734 * old - 0.0022))**(-1)
          xnew = old - funct/fprime
          if (abs(xnew - old).lt.tol) return
          old = xnew
        end do

      endif

      print *,'max iteration exceeded when solving for Utau, Utau=',old

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getsublim(z,rh,tair,Utau,
     &  z_0,V_susp,V_salt,Utau_t,ht_rhobs,flag,pi)

      implicit none

      real pi,ro_ice,xM,R,R_dryair,vonKarman,visc_air,h_s
      real xlamdaT,D,ro_sat,rh_offset,sigma
      real alfa,rbar_r,xmbar,rbar,u_z,x_r,wbar,flag
      real z,rh,tair,Utau,z_0,V_susp,V_salt,Utau_t,ht_rhobs
      real V_r,xN_r,xNu,xSh,tmp1,tmp2,top,bottom,V_rsalt

      ro_ice = 917.0
      xM = 18.01
      R = 8313.
      R_dryair = 287.
      vonKarman = 0.4
      visc_air = 13.e-6
      h_s = 2.838e6

c     xlamdaT = 0.00063 * tair + 0.0673
      xlamdaT = 0.024
      D = 2.06e-5 * (tair/273.)**(1.75)
c     ro_sat = 0.622 * 10.0**(11.40 - 2353./tair) / (R_dryair * tair)
      ro_sat = 0.622 / (R_dryair * tair) *
     &  610.78 * exp(21.875 * (tair - 273.16) / (tair - 7.66))

c Assume that the rh varies according to a modification to 
c   Pomeroy's humidity variation with height equation.
      rh_offset = 1.0 - 0.027 * log(ht_rhobs)
      sigma = (0.01 * rh - 1.0) * (rh_offset + 0.027 * log(z))
      sigma = min(0.0,sigma)
      sigma = max(-1.0,sigma)

      alfa = 4.08 + 12.6 * z
      rbar_r = 4.6e-5 * z**(-0.258)
      xmbar = 4.0/3.0 * pi * ro_ice * rbar_r**3 *
     &  (1.0 + 3.0/alfa + 2.0/alfa**2)
      rbar = ((3.0 * xmbar) / (4.0 * pi * ro_ice))**(0.33)
      u_z = Utau/vonKarman * log(z/z_0)
      x_r = 0.005 * u_z**(1.36)
      wbar = 1.1e7 * rbar**(1.8)

      if (flag.eq.1.0) then

c Compute the sublimation loss rate coefficient for the suspension
c   layer.
        V_r = wbar + 3.0 * x_r * cos(pi/4.0)
        xN_r = 2.0 * rbar * V_r / visc_air
        xNu = 1.79 + 0.606 * xN_r**(0.5)
        xSh = xNu
        tmp1 = (h_s * xM)/(R * tair) - 1.0
        tmp2 = xlamdaT * tair * xNu
        top = 2.0 * pi * rbar * sigma
        bottom = h_s/tmp2 * tmp1 + 1.0/(D * ro_sat * xSh)
        V_susp = (top/bottom)/xmbar
        V_salt = 0.0

      elseif (flag.eq.0.0) then

c Compute the sublimation loss rate coefficient for the saltation
c   layer.
        V_rsalt = 0.68 * Utau + 2.3 * Utau_t
        xN_r = 2.0 * rbar * V_rsalt / visc_air
        xNu = 1.79 + 0.606 * xN_r**(0.5)
        xSh = xNu
        tmp1 = (h_s * xM)/(R * tair) - 1.0
        tmp2 = xlamdaT * tair * xNu
        top = 2.0 * pi * rbar * sigma
        bottom = h_s/tmp2 * tmp1 + 1.0/(D * ro_sat * xSh)
        V_salt = (top/bottom)/xmbar
        V_susp = 0.0

      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getdirection(nx,ny,uwind_grid,vwind_grid,index_ue,
     &  index_uw,index_vn,index_vs)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,npairs

      real sign1,sign2

      integer index_ue(ny_max,2*nx_max+1)
      integer index_uw(ny_max,2*nx_max+1)
      integer index_vn(nx_max,2*ny_max+1)
      integer index_vs(nx_max,2*ny_max+1)

      real uwind_grid(nx_max,ny_max)
      real vwind_grid(nx_max,ny_max)

c Index whether the winds are blowing east or west.  The first
c   column of the index array is the number of pairs of begining
c   and ending array index of blocks of wind running in the same
c   direction.

c Sweep looking for WESTERLY winds, looking for positive numbers.
      do j=1,ny

        if (uwind_grid(1,j).le.0.0) then
          sign1 = -1.0
        else
          sign1 = 1.0
        endif

          if (sign1.gt.0.0) then
            npairs = 1
            index_uw(j,2) =  1
          else
            npairs = 0
          endif
        do i=2,nx

          if (uwind_grid(i-1,j).le.0.0) then
            sign1 = -1.0
          else
            sign1 = 1.0
          endif
          if (uwind_grid(i,j).le.0.0) then
            sign2 = -1.0
          else
            sign2 = 1.0
          endif

          if (sign2.ne.sign1) then
c We have a sign change.
            if (sign2.gt.0.0) then
c We have gone from negative to positive, indicating the start
c   of a new positive group.
              npairs = npairs + 1
              index_uw(j,npairs*2) = i
            else
c We have gone from positive to negative, indicating the end of
c   the group.
              index_uw(j,npairs*2+1) = i - 1
            endif
          endif
        enddo

        if (uwind_grid(nx,j).le.0.0) then
          sign1 = -1.0
        else
          sign1 = 1.0
        endif

        if (sign1.gt.0.0) then
          index_uw(j,npairs*2+1) = nx
        endif
        index_uw(j,1) = npairs
      enddo

c     do j=1,ny
c       print 30, (index_uw(j,k),k=1,index_uw(j,1)*2+1)
c     enddo
c     print *
c     print *

c Sweep looking for EASTERLY winds, looking for negative numbers.
      do j=1,ny

        if (uwind_grid(1,j).le.0.0) then
          sign1 = -1.0
        else
          sign1 = 1.0
        endif

          if (sign1.lt.0.0) then
            npairs = 1
            index_ue(j,2) = 1
          else
            npairs = 0
          endif
        do i=2,nx

          if (uwind_grid(i-1,j).le.0.0) then
            sign1 = -1.0
          else
            sign1 = 1.0
          endif
          if (uwind_grid(i,j).le.0.0) then
            sign2 = -1.0
          else
            sign2 = 1.0
          endif

          if (sign2.ne.sign1) then
c We have a sign change.
            if (sign2.lt.0.0) then
c We have gone from positive to negative, indicating the start
c   of a new negative group.
              npairs = npairs + 1
              index_ue(j,npairs*2) = i
            else
c We have gone from negative to positive, indicating the end of
c   the group.
              index_ue(j,npairs*2+1) = i - 1
            endif
          endif
        enddo

        if (uwind_grid(nx,j).le.0.0) then
          sign1 = -1.0
        else
          sign1 = 1.0
        endif

        if (sign1.lt.0.0) then
          index_ue(j,npairs*2+1) = nx
        endif
        index_ue(j,1) = npairs
      enddo

c     do j=1,ny
c       print 30, (index_ue(j,k),k=1,index_ue(j,1)*2+1)
c     enddo
c     print *
c     print *

c Sweep looking for SOUTHERLY winds, looking for positive numbers.
      do i=1,nx

        if (vwind_grid(i,1).le.0.0) then
          sign1 = -1.0
        else
          sign1 = 1.0
        endif

          if (sign1.gt.0.0) then
            npairs = 1
            index_vs(i,2) = 1
          else
            npairs = 0
          endif
        do j=2,ny

          if (vwind_grid(i,j-1).le.0.0) then
            sign1 = -1.0
          else
            sign1 = 1.0
          endif
          if (vwind_grid(i,j).le.0.0) then
            sign2 = -1.0
          else
            sign2 = 1.0
          endif

          if (sign2.ne.sign1) then
c We have a sign change.
            if (sign2.gt.0.0) then
c We have gone from negative to positive, indicating the start
c   of a new positive group.
              npairs = npairs + 1
              index_vs(i,npairs*2) = j
            else
c We have gone from positive to negative, indicating the end of
c   the group.
              index_vs(i,npairs*2+1) = j - 1
            endif
          endif
        enddo

        if (vwind_grid(i,ny).le.0.0) then
          sign1 = -1.0
        else
          sign1 = 1.0
        endif

        if (sign1.gt.0.0) then
          index_vs(i,npairs*2+1) = ny
        endif
        index_vs(i,1) = npairs
      enddo

c     do i=1,nx
c       print 30, (index_vs(i,k),k=1,index_vs(i,1)*2+1)
c     enddo
c     print *
c     print *

c Sweep looking for NORTHERLY winds, looking for negative numbers.
      do i=1,nx

        if (vwind_grid(i,1).le.0.0) then
          sign1 = -1.0
        else
          sign1 = 1.0
        endif

          if (sign1.lt.0.0) then
            npairs = 1
            index_vn(i,2) = 1
          else
            npairs = 0
          endif
        do j=2,ny

          if (vwind_grid(i,j-1).le.0.0) then
            sign1 = -1.0
          else
            sign1 = 1.0
          endif
          if (vwind_grid(i,j).le.0.0) then
            sign2 = -1.0
          else
            sign2 = 1.0
          endif

          if (sign2.ne.sign1) then
c We have a sign change.
            if (sign2.lt.0.0) then
c We have gone from positive to negative, indicating the start
c   of a new negative group.
              npairs = npairs + 1
              index_vn(i,npairs*2) = j
            else
c We have gone from negative to positive, indicating the end of
c   the group.
              index_vn(i,npairs*2+1) = j - 1
            endif
          endif
        enddo

        if (vwind_grid(i,ny).le.0.0) then
          sign1 = -1.0
        else
          sign1 = 1.0
        endif

        if (sign1.lt.0.0) then
          index_vn(i,npairs*2+1) = ny
        endif
        index_vn(i,1) = npairs
      enddo

c     do i=1,nx
c       print 30, (index_vn(i,k),k=1,index_vn(i,1)*2+1)
c     enddo
c     print *
c     print *
c 30  format(20i4)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getnewdepth(nx,ny,deltax,deltay,Qsalt_u,
     &  Qsalt_v,dh_salt_u,dh_salt_v,index_ue,index_uw,
     &  index_vn,index_vs,ro_snow,dt,vegsnowd_xy,snow_d,
     &  soft_snow_d)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny
      integer k,istart,iend,jstart,jend

      real deltax,deltay,ro_snow,dt,dQsalt,snowdmin
      real hard_snow_d,weight_u,weight_v,eps

      integer index_ue(ny_max,2*nx_max+1)
      integer index_uw(ny_max,2*nx_max+1)
      integer index_vn(nx_max,2*ny_max+1)
      integer index_vs(nx_max,2*ny_max+1)

      real Qsalt_u(nx_max,ny_max)
      real Qsalt_v(nx_max,ny_max)

      real snow_d(nx_max,ny_max)
      real soft_snow_d(nx_max,ny_max)
      real dh_salt_u(nx_max,ny_max)
      real dh_salt_v(nx_max,ny_max)
      real vegsnowd_xy(nx_max,ny_max)

c Define an upwind boundary condition for saltation (here I have
c   assumed that the transport is in equilibrium).
      do i=1,nx
        do j=1,ny
          dh_salt_u(i,j) = 0.0
          dh_salt_v(i,j) = 0.0
        enddo
      enddo

c Consider WESTERLY winds.
      do j=1,ny
        do k=1,index_uw(j,1)
          istart = index_uw(j,k*2)+1
          iend = index_uw(j,k*2+1)
          do i=istart,iend
            dQsalt = Qsalt_u(i,j) - Qsalt_u(i-1,j)
            dh_salt_u(i,j) = (- dt) / ro_snow * dQsalt / deltax

c Make adjustments for the case where there is no snow available
c   on the ground (or captured within the vegetation) to be
c   eroded.
            hard_snow_d = snow_d(i,j) - soft_snow_d(i,j)
            snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)
            if (snow_d(i,j).gt.snowdmin) then
              if (snow_d(i,j)+dh_salt_u(i,j).le.snowdmin) then
                dh_salt_u(i,j) = snowdmin - snow_d(i,j)
                Qsalt_u(i,j) = Qsalt_u(i-1,j) - dh_salt_u(i,j) *
     &            ro_snow * deltax / dt
              endif
            else
              Qsalt_u(i,j) = 0.0
              dh_salt_u(i,j) = 0.0
            endif
          enddo
        enddo
      enddo

c Consider EASTERLY winds.
      do j=1,ny
        do k=1,index_ue(j,1)
          iend = index_ue(j,k*2)
          istart = index_ue(j,k*2+1)-1
          do i=istart,iend,-1
            dQsalt = Qsalt_u(i,j) - Qsalt_u(i+1,j)
            dh_salt_u(i,j) = (- dt) / ro_snow * dQsalt / deltax

c Make adjustments for the case where there is no snow available
c   on the ground (or captured within the vegetation) to be
c   eroded.
            hard_snow_d = snow_d(i,j) - soft_snow_d(i,j)
            snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)
            if (snow_d(i,j).gt.snowdmin) then
              if (snow_d(i,j)+dh_salt_u(i,j).le.snowdmin) then
                dh_salt_u(i,j) = snowdmin - snow_d(i,j)
                Qsalt_u(i,j) = Qsalt_u(i+1,j) - dh_salt_u(i,j) *
     &            ro_snow * deltax / dt
              endif
            else
              Qsalt_u(i,j) = 0.0
              dh_salt_u(i,j) = 0.0
            endif
          enddo
        enddo
      enddo

c Consider SOUTHERLY winds.
      do i=1,nx
        do k=1,index_vs(i,1)
          jstart = index_vs(i,k*2)+1
          jend = index_vs(i,k*2+1)
          do j=jstart,jend
            dQsalt = Qsalt_v(i,j) - Qsalt_v(i,j-1)
            dh_salt_v(i,j) = (- dt) / ro_snow * dQsalt / deltay

c Make adjustments for the case where there is no snow available
c   on the ground (or captured within the vegetation) to be
c   eroded.
            hard_snow_d = snow_d(i,j) - soft_snow_d(i,j)
            snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)
            if (snow_d(i,j).gt.snowdmin) then
              if (snow_d(i,j)+dh_salt_v(i,j).le.snowdmin) then
                dh_salt_v(i,j) = snowdmin - snow_d(i,j)
                Qsalt_v(i,j) = Qsalt_v(i,j-1) - dh_salt_v(i,j) *
     &            ro_snow * deltay / dt
              endif
            else
              Qsalt_v(i,j) = 0.0
              dh_salt_v(i,j) = 0.0
            endif
          enddo
        enddo
      enddo

c Consider NORTHERLY winds.
      do i=1,nx
        do k=1,index_vn(i,1)
          jend = index_vn(i,k*2)
          jstart = index_vn(i,k*2+1)-1
          do j=jstart,jend,-1
            dQsalt = Qsalt_v(i,j) - Qsalt_v(i,j+1)
            dh_salt_v(i,j) = (- dt) / ro_snow * dQsalt / deltay

c Make adjustments for the case where there is no snow available
c   on the ground (or captured within the vegetation) to be
c   eroded.
            hard_snow_d = snow_d(i,j) - soft_snow_d(i,j)
            snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)
            if (snow_d(i,j).gt.snowdmin) then
              if (snow_d(i,j)+dh_salt_v(i,j).le.snowdmin) then
                dh_salt_v(i,j) = snowdmin - snow_d(i,j)
                Qsalt_v(i,j) = Qsalt_v(i,j+1) - dh_salt_v(i,j) *
     &            ro_snow * deltay / dt
              endif
            else
              Qsalt_v(i,j) = 0.0
              dh_salt_v(i,j) = 0.0
            endif
          enddo
        enddo
      enddo

c Update the snow depth changes due to saltation transport from the
c   the east and west, and north and south.  Also correct dh_salt_u
c   and dh_salt_v to account for the minimum snow depth.
      eps = 1e-6
      do i=1,nx
        do j=1,ny
          weight_u = abs(dh_salt_u(i,j)) /
     &      (abs(dh_salt_u(i,j)) + abs(dh_salt_v(i,j)) + eps)

          weight_v = abs(dh_salt_v(i,j)) /
     &      (abs(dh_salt_u(i,j)) + abs(dh_salt_v(i,j)) + eps)

          dh_salt_u(i,j) = weight_u * dh_salt_u(i,j)
          dh_salt_v(i,j) = weight_v * dh_salt_v(i,j)

          snow_d(i,j) = snow_d(i,j) + dh_salt_u(i,j) + dh_salt_v(i,j)

          soft_snow_d(i,j) = soft_snow_d(i,j) + dh_salt_u(i,j) +
     &      dh_salt_v(i,j)
        enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine subgrid_1(nx,ny,snow_d,
     &  index_ue,index_uw,index_vn,index_vs,
     &  tabler_nn,tabler_ss,tabler_ee,tabler_ww,
     &  tabler_ne,tabler_se,tabler_sw,tabler_nw,uwind_grid,
     &  vwind_grid,tabler_dir)

c This subroutine forces SnowTran-3D's snow accumluation profiles
c   to be bounded by the equilibrium topographic drift catchment
c   profiles observed and modeled by Tabler (1975).

c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
c   topographic catchments.  Proceedings of the 43rd Annual Western
c   Snow Conference, San Diego, California, 87-97.

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny
      integer k,istart,iend,jstart,jend

      real snow_d_extra,snow_sfc,tabler,tabler_dir

      integer index_ue(ny_max,2*nx_max+1)
      integer index_uw(ny_max,2*nx_max+1)
      integer index_vn(nx_max,2*ny_max+1)
      integer index_vs(nx_max,2*ny_max+1)

      real snow_d(nx_max,ny_max)
      real snow_d1(nx_max,ny_max)
      real snow_d2(nx_max,ny_max)
      real tabler_nn(nx_max,ny_max)
      real tabler_ss(nx_max,ny_max)
      real tabler_ee(nx_max,ny_max)
      real tabler_ww(nx_max,ny_max)
      real tabler_ne(nx_max,ny_max)
      real tabler_se(nx_max,ny_max)
      real tabler_sw(nx_max,ny_max)
      real tabler_nw(nx_max,ny_max)
      real uwind_grid(nx_max,ny_max)
      real vwind_grid(nx_max,ny_max)

      real weight_u(nx_max,ny_max)
      real weight_v(nx_max,ny_max)

c This is just a summary of all of the possibilities.
c           if(winddir(i,j).gt.337.5.and.winddir(i,j).le.22.5)then
c             tabler = tabler_nn(i,j)
c           elseif(winddir(i,j).gt.22.5.and.winddir(i,j).le.67.5)then
c             tabler = tabler_ne(i,j)
c           elseif(winddir(i,j).gt.67.5.and.winddir(i,j).le.112.5)then
c             tabler = tabler_ee(i,j)
c           elseif(winddir(i,j).gt.112.5.and.winddir(i,j).le.157.5)then
c             tabler = tabler_se(i,j)
c           elseif(winddir(i,j).gt.157.5.and.winddir(i,j).le.202.5)then
c             tabler = tabler_ss(i,j)
c           elseif(winddir(i,j).gt.202.5.and.winddir(i,j).le.247.5)then
c             tabler = tabler_sw(i,j)
c           elseif(winddir(i,j).gt.247.5.and.winddir(i,j).le.292.5)then
c             tabler = tabler_ww(i,j)
c           elseif(winddir(i,j).gt.292.5.and.winddir(i,j).le.337.5)then
c             tabler = tabler_nw(i,j)
c           endif

c Create a copy of the incoming snow depth distribution.  Also define
c   the u and v weighting functions.
      do j=1,ny
        do i=1,nx
          snow_d1(i,j) = snow_d(i,j)
          snow_d2(i,j) = snow_d(i,j)

          weight_u(i,j) = abs(uwind_grid(i,j)) /
     &          sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
          weight_v(i,j) = abs(vwind_grid(i,j)) /
     &          sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
        enddo
      enddo

c Consider WESTERLY winds.
      do j=1,ny
        do k=1,index_uw(j,1)
          istart = index_uw(j,k*2)+1
          iend = index_uw(j,k*2+1)
          do i=istart,iend

            if(tabler_dir.gt.337.5.and.tabler_dir.le.360.0.or.
     &        tabler_dir.ge.0.0.and.tabler_dir.le.22.5)then
              tabler = tabler_nn(i,j)
            elseif(tabler_dir.gt.157.5.and.tabler_dir.le.202.5)then
              tabler = tabler_ss(i,j)
            elseif(tabler_dir.gt.202.5.and.tabler_dir.le.247.5)then
              tabler = tabler_sw(i,j)
            elseif(tabler_dir.gt.247.5.and.tabler_dir.le.292.5)then
              tabler = tabler_ww(i,j)
            elseif(tabler_dir.gt.292.5.and.tabler_dir.le.337.5)then
              tabler = tabler_nw(i,j)
            endif

            snow_sfc = tabler

            if (snow_d1(i,j).gt.snow_sfc) then
              snow_d_extra = (snow_d1(i,j) - snow_sfc) * weight_u(i,j)
              snow_d1(i,j) = snow_d1(i,j) - snow_d_extra
              if (i.lt.nx) then
                snow_d1(i+1,j) = snow_d1(i+1,j) + snow_d_extra
              else
                snow_d1(i,j) = snow_d1(i,j)
              endif
            endif

          enddo
        enddo
      enddo

c Consider EASTERLY winds.
      do j=1,ny
        do k=1,index_ue(j,1)
          iend = index_ue(j,k*2)
          istart = index_ue(j,k*2+1)-1
          do i=istart,iend,-1

            if(tabler_dir.gt.337.5.and.tabler_dir.le.360.0.or.
     &        tabler_dir.ge.0.0.and.tabler_dir.le.22.5)then
              tabler = tabler_nn(i,j)
            elseif(tabler_dir.gt.22.5.and.tabler_dir.le.67.5)then
              tabler = tabler_ne(i,j)
            elseif(tabler_dir.gt.67.5.and.tabler_dir.le.112.5)then
              tabler = tabler_ee(i,j)
            elseif(tabler_dir.gt.112.5.and.tabler_dir.le.157.5)then
              tabler = tabler_se(i,j)
            elseif(tabler_dir.gt.157.5.and.tabler_dir.le.202.5)then
              tabler = tabler_ss(i,j)
            endif

            snow_sfc = tabler

            if (snow_d1(i,j).gt.snow_sfc) then
              snow_d_extra = (snow_d1(i,j) - snow_sfc) * weight_u(i,j)
              snow_d1(i,j) = snow_d1(i,j) - snow_d_extra
              if (i.gt.1) then
                snow_d1(i-1,j) = snow_d1(i-1,j) + snow_d_extra
              else
                snow_d1(i,j) = snow_d1(i,j)
              endif
            endif
          enddo
        enddo
      enddo

c Consider SOUTHERLY winds.
      do i=1,nx
        do k=1,index_vs(i,1)
          jstart = index_vs(i,k*2)+1
          jend = index_vs(i,k*2+1)
          do j=jstart,jend

            if(tabler_dir.gt.67.5.and.tabler_dir.le.112.5)then
              tabler = tabler_ee(i,j)
            elseif(tabler_dir.gt.112.5.and.tabler_dir.le.157.5)then
              tabler = tabler_se(i,j)
            elseif(tabler_dir.gt.157.5.and.tabler_dir.le.202.5)then
              tabler = tabler_ss(i,j)
            elseif(tabler_dir.gt.202.5.and.tabler_dir.le.247.5)then
              tabler = tabler_sw(i,j)
            elseif(tabler_dir.gt.247.5.and.tabler_dir.le.292.5)then
              tabler = tabler_ww(i,j)
            endif

            snow_sfc = tabler

            if (snow_d2(i,j).gt.snow_sfc) then
              snow_d_extra = (snow_d2(i,j) - snow_sfc) * weight_v(i,j)
              snow_d2(i,j) = snow_d2(i,j) - snow_d_extra
              if (j.lt.ny) then
                snow_d2(i,j+1) = snow_d2(i,j+1) + snow_d_extra
              else
                snow_d2(i,j) = snow_d2(i,j)
              endif
            endif
          enddo
        enddo
      enddo

c Consider NORTHERLY winds.
      do i=1,nx
        do k=1,index_vn(i,1)
          jend = index_vn(i,k*2)
          jstart = index_vn(i,k*2+1)-1
          do j=jstart,jend,-1

            if(tabler_dir.gt.337.5.and.tabler_dir.le.360.0.or.
     &        tabler_dir.ge.0.0.and.tabler_dir.le.22.5)then
              tabler = tabler_nn(i,j)
            elseif(tabler_dir.gt.22.5.and.tabler_dir.le.67.5)then
              tabler = tabler_ne(i,j)
            elseif(tabler_dir.gt.67.5.and.tabler_dir.le.112.5)then
              tabler = tabler_ee(i,j)
            elseif(tabler_dir.gt.247.5.and.tabler_dir.le.292.5)then
              tabler = tabler_ww(i,j)
            elseif(tabler_dir.gt.292.5.and.tabler_dir.le.337.5)then
              tabler = tabler_nw(i,j)
            endif

            snow_sfc = tabler

            if (snow_d2(i,j).gt.snow_sfc) then
              snow_d_extra = (snow_d2(i,j) - snow_sfc) * weight_v(i,j)
              snow_d2(i,j) = snow_d2(i,j) - snow_d_extra
              if (j.gt.1) then
                snow_d2(i,j-1) = snow_d2(i,j-1) + snow_d_extra
              else
                snow_d2(i,j) = snow_d2(i,j)
              endif
            endif

          enddo
        enddo
      enddo

c Update the snow depths resulting from these redistributions.
      do j=1,ny
        do i=1,nx
          snow_d(i,j) = snow_d1(i,j) * weight_u(i,j) +
     &      snow_d2(i,j) * weight_v(i,j)
        enddo
      enddo

c Clean up the boundaries.  Make the boundary values equal to
c   the values just inside the boundaries.
      do i=2,nx-1
        snow_d(i,1) = snow_d(i,2)
        snow_d(i,ny) = snow_d(i,ny-1)
      enddo
      do j=1,ny
        snow_d(1,j) = snow_d(2,j)
        snow_d(nx,j) = snow_d(nx-1,j)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tabler_3d(nx,ny,topo_land,deltax,deltay,
     &  tabler_ww,tabler_ee,tabler_ss,tabler_nn,erosion_dist,
     &  tabler_ne,tabler_se,tabler_sw,tabler_nw,slope_adjust)

c This subroutine uses Tabler (1975) to define equilibrium profiles
c   for the topographic drift catchments, for the case of winds
c   from the EAST, WEST, NORTH, and SOUTH, and anywhere inbetween.

      implicit none

      include 'snowmodel.inc'

      integer nx,ny,irotate_flag

      real deltax,deltay,erosion_dist,slope_adjust

      real topo_land(nx_max,ny_max)

      real tabler_nn(nx_max,ny_max)
      real tabler_ss(nx_max,ny_max)
      real tabler_ee(nx_max,ny_max)
      real tabler_ww(nx_max,ny_max)

      real tabler_ne(nx_max,ny_max)
      real tabler_se(nx_max,ny_max)
      real tabler_sw(nx_max,ny_max)
      real tabler_nw(nx_max,ny_max)

c Here we generate maximum snow accumulation surfaces for n, ne, e,
c   se, s, sw, w, and nw winds.  I call these "tabler surfaces".  
c
c They are valid for the wind direction ranges: N=337.5-22.5,
c   NE=22.5-67.5, E=67.5-112.5, SE=112.5-157.5, S=157.5-202.5,
c   SW=202.5-247.5, W=247.5-292.5, and NW=292.5-337.5.
c
c These Tabler Surfaces define a "potential" snow surface that
c   represents the maximum possible snow-accumulation depth from winds
c   coming from these directions.
c
c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
c   topographic catchments.  Proceedings of the 43rd Annual Western
c   Snow Conference, San Diego, California, 87-97.

c Consider N winds.
      irotate_flag = 1
      call tabler_n(nx,ny,topo_land,tabler_nn,deltay,
     &  erosion_dist,irotate_flag,slope_adjust)

c Consider NE winds.
      irotate_flag = 2
      call tabler_e(nx,ny,topo_land,tabler_ne,1.41*deltax,
     &  erosion_dist,irotate_flag,slope_adjust)

c Consider E winds.
      irotate_flag = 1
      call tabler_e(nx,ny,topo_land,tabler_ee,deltax,
     &  erosion_dist,irotate_flag,slope_adjust)

c Consider SE winds.
      irotate_flag = 2
      call tabler_s(nx,ny,topo_land,tabler_se,1.41*deltay,
     &  erosion_dist,irotate_flag,slope_adjust)

c Consider S winds.
      irotate_flag = 1
      call tabler_s(nx,ny,topo_land,tabler_ss,deltay,
     &  erosion_dist,irotate_flag,slope_adjust)

c Consider SW winds.
      irotate_flag = 2
      call tabler_w(nx,ny,topo_land,tabler_sw,1.41*deltax,
     &  erosion_dist,irotate_flag,slope_adjust)

c Consider W winds.
      irotate_flag = 1
      call tabler_w(nx,ny,topo_land,tabler_ww,deltax,
     &  erosion_dist,irotate_flag,slope_adjust)

c Consider NW winds.
      irotate_flag = 2
      call tabler_n(nx,ny,topo_land,tabler_nw,1.41*deltay,
     &  erosion_dist,irotate_flag,slope_adjust)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tabler_w(nx,ny,topo_land,tabler_ww,deltax,
     &  erosion_dist,irotate_flag,slope_adjust)

c This subroutine uses Tabler (1975) to define equilibrium profiles
c   for the topographic drift catchments, for the case of winds
c   from the WEST and SOUTHWEST.

c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
c   topographic catchments.  Proceedings of the 43rd Annual Western
c   Snow Conference, San Diego, California, 87-97.

      implicit none

      include 'snowmodel.inc'

      integer nx,ny,i,j,istart,iend,JJ,irotate_flag,nny,nnx,
     &  ii,iii,nx_max_tabler,maxlines
      integer npts
      real deltax,y,erosion_dist,slope_adjust,xmax_slope,dx,test,
     &  x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

      real topo_land(nx_max,ny_max)
      real tabler_ww(nx_max,ny_max)
      real tabler(nx_max)
      real topo_line(nx_max)
      real drift_start_topo(nx_max)

      parameter (nx_max_tabler=nx_max*1000)
      real topo_1m(nx_max_tabler)
      real tabler_1m(nx_max_tabler)
      real drift_start_topo_1m(nx_max_tabler)

c This program:
c   1) takes the coarse model topography and generates 1.0-m grid
c        increment topo lines;
c   2) uses those to generate the Tabler surface profiles; and
c   3) extracts the profiles at the coarse model grid cells.

c Fill the snow array with topo data.
      do j=1,ny
        do i=1,nx
          tabler_ww(i,j) = topo_land(i,j)
        enddo
      enddo

c Define the length of the j-looping, depending on whether there
c   is any rotation done to get 45 deg winds.
      if (irotate_flag.eq.2) then
        nny = nx+ny-1
      else
        nny = ny
      endif

c Required parameters.
      if (irotate_flag.eq.2) then
        dx = 1.41
        test = amod(deltax/1.41,dx/1.41)
      else
        dx = 1.0
        test = amod(deltax,dx)
      endif

      xmax_slope = -0.20

c This Tabler program has not been made general enough to deal
c   with deltax and deltay values that are not evenly divisible
c   by 1.0 (or 1.41 for diagonal profiles).
      if (abs(test).gt.1.0e10-5) then
        print *,'To generate the Tabler surfaces, deltax and deltay'
        print *,'  must be evenly divisible by 1.0 or 1.41.'
        print *,'  deltax = ',deltax
        stop
      endif

c Define the number of 1-m grid cells in each model grid cell, and
c   calculate how many of these are in the entire nx domain.
      nnx = nint(deltax/dx)
      maxlines = (nx - 1) * nnx + 1

c Define the starting and ending points of the line we are going to
c   work with.  This is done to make sure we are never looking
c   outside the data for numbers to work with.
      istart = 1
      if (irotate_flag.eq.2) then
        iend = maxlines - (32+11+10+11)
      else
        iend = maxlines - (45+15+15+15)
      endif

c Extract the line we are going to work with.
      do j=1,nny

        if (irotate_flag.eq.2) then
          do i=1,nx
            JJ = j + i - nx
            drift_start_topo(i) = 0.0
            if (JJ.le.0) then
              tabler(i) = tabler_ww(1-j+nx,1)
              topo_line(i) = tabler_ww(1-j+nx,1)
            elseif (JJ.gt.ny) then
              tabler(i) = tabler(i-1)
              topo_line(i) = tabler(i-1)
            else
              tabler(i) = tabler_ww(i,JJ)
              topo_line(i) = tabler_ww(i,JJ)
            endif
          enddo
        else
          do i=1,nx
            tabler(i) = tabler_ww(i,j)
            topo_line(i) = topo_land(i,j)
            drift_start_topo(i) = 0.0
          enddo
        endif

c To build the 1.0 m line, use linear interpolation between the
c   model topo data.  Include the end point.
        do i=1,nx-1
          do ii=1,nnx
            iii = (i - 1) * nnx + ii
            x1 = 0.0
            x = real(ii - 1) * dx
            y2 = topo_line(i+1)
            y1 = topo_line(i)
            topo_1m(iii) = y1 + ((y2 - y1)/deltax) * (x - x1)
          enddo
        enddo
        topo_1m((nx - 1) * nnx + 1) = topo_line(nx)

c Use this topo array to be the starting point for generating the
c   Tabler surfaces.
        do i=1,maxlines
          tabler_1m(i) = topo_1m(i)
          drift_start_topo_1m(i) = 0.0
        enddo

c Run the Tabler model.
        do i=istart,iend
          if (irotate_flag.eq.2) then
            t1 = tabler_1m(i)
            t2 = tabler_1m(i+31)
            t3 = tabler_1m(i+31+11)
            t4 = tabler_1m(i+31+21)
            t5 = tabler_1m(i+31+32)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(i+32) = max(topo_1m(i+32),
     &        tabler_1m(i+31) + y * slope_adjust * dx)

c Create an erosion area 'erosion_dist' upwind of the tabler surface.
            npts = nint(erosion_dist/1.41)
            if (tabler_1m(i+32).ne.topo_1m(i+32)) then
              do ii=1,npts
                if (tabler_1m(i+32-ii).eq.topo_1m(i+32-ii))
     &            drift_start_topo_1m(i+32-ii) = -8888.0
              enddo
            endif
          else
            t1 = tabler_1m(i)
            t2 = tabler_1m(i+44)
            t3 = tabler_1m(i+44+15)
            t4 = tabler_1m(i+44+30)
            t5 = tabler_1m(i+44+45)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(i+45) = max(topo_1m(i+45),
     &        tabler_1m(i+44) + y * slope_adjust * dx)

c Create an erosion area 'erosion_dist' upwind of the tabler surface.
            npts = nint(erosion_dist/1.0)
            if (tabler_1m(i+45).ne.topo_1m(i+45)) then
              do ii=1,npts
                if (tabler_1m(i+45-ii).eq.topo_1m(i+45-ii))
     &            drift_start_topo_1m(i+45-ii) = -8888.0
              enddo
            endif
          endif
        enddo

c Extract the profile at the model grid points.
        do i=1,nx
          ii = (i - 1) * nnx + 1
          tabler(i) = tabler_1m(ii)
          drift_start_topo(i) = drift_start_topo_1m(ii)
        enddo

c Use the 1-D arrays to fill in the 2-D tabler-surface array.
        do i=1,nx
          if (irotate_flag.eq.2) then
            JJ = j + i - nx
            if (JJ.ge.1 .and. JJ.le.ny) then
              tabler_ww(i,JJ) = tabler(i) + drift_start_topo(i)
            endif
          else
            tabler_ww(i,j) = tabler(i) + drift_start_topo(i)
          endif
        enddo

      enddo

c Convert snow_traps back to actual snow depths instead of
c   depth plus topography.
      do j=1,ny
        do i=1,nx
          tabler_ww(i,j) = tabler_ww(i,j) - topo_land(i,j)
        enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tabler_e(nx,ny,topo_land,tabler_ee,deltax,
     &  erosion_dist,irotate_flag,slope_adjust)

c This subroutine uses Tabler (1975) to define equilibrium profiles
c   for the topographic drift catchments, for the case of winds
c   from the EAST and NORTHEAST.

c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
c   topographic catchments.  Proceedings of the 43rd Annual Western
c   Snow Conference, San Diego, California, 87-97.

      implicit none

      include 'snowmodel.inc'

      integer nx,ny,i,j,istart,iend,JJ,irotate_flag,nny,nnx,
     &  ii,iii,nx_max_tabler,maxlines
      integer npts
      real deltax,y,erosion_dist,slope_adjust,xmax_slope,dx,test,
     &  x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

      real topo_land(nx_max,ny_max)
      real tabler_ee(nx_max,ny_max)
      real tabler(nx_max)
      real topo_line(nx_max)
      real drift_start_topo(nx_max)

      parameter (nx_max_tabler=nx_max*1000)
      real topo_1m(nx_max_tabler)
      real tabler_1m(nx_max_tabler)
      real drift_start_topo_1m(nx_max_tabler)

c This program:
c   1) takes the coarse model topography and generates 1.0-m grid
c        increment topo lines;
c   2) uses those to generate the Tabler surface profiles; and
c   3) extracts the profiles at the coarse model grid cells.

c Fill the snow array with topo data.
      do j=1,ny
        do i=1,nx
          tabler_ee(i,j) = topo_land(i,j)
        enddo
      enddo

c Define the length of the j-looping, depending on whether there
c   is any rotation done to get 45 deg winds.
      if (irotate_flag.eq.2) then
        nny = nx+ny-1
      else
        nny = ny
      endif

c Required parameters.
      if (irotate_flag.eq.2) then
        dx = 1.41
        test = amod(deltax/1.41,dx/1.41)
      else
        dx = 1.0
        test = amod(deltax,dx)
      endif

      xmax_slope = -0.20

c This Tabler program has not been made general enough to deal
c   with deltax and deltay values that are not evenly divisible
c   by 1.0 (or 1.41 for diagonal profiles).
      if (abs(test).gt.1.0e10-5) then
        print *,'To generate the Tabler surfaces, deltax and deltay'
        print *,'  must be evenly divisible by 1.0 or 1.41.'
        print *,'  deltax = ',deltax
        stop
      endif

c Define the number of 1-m grid cells in each model grid cell, and
c   calculate how many of these are in the entire nx domain.
      nnx = nint(deltax/dx)
      maxlines = (nx - 1) * nnx + 1

c Define the starting and ending points of the line we are going to
c   work with.  This is done to make sure we are never looking
c   outside the data for numbers to work with.
      istart = maxlines
      if (irotate_flag.eq.2) then
        iend = 1 + (32+11+10+11)
      else
        iend = 1 + (45+15+15+15)
      endif

c Extract the line we are going to work with.
      do j=1,nny

        if (irotate_flag.eq.2) then
          do i=1,nx
            JJ = j + i - nx
            drift_start_topo(i) = 0.0
            if (JJ.le.0) then
              tabler(i) = tabler_ee(1-j+nx,1)
              topo_line(i) = tabler_ee(1-j+nx,1)
            elseif (JJ.gt.ny) then
              tabler(i) = tabler(i-1)
              topo_line(i) = tabler(i-1)
            else
              tabler(i) = tabler_ee(i,JJ)
              topo_line(i) = tabler_ee(i,JJ)
            endif
          enddo
        else
          do i=1,nx
            tabler(i) = tabler_ee(i,j)
            topo_line(i) = topo_land(i,j)
            drift_start_topo(i) = 0.0
          enddo
        endif

c To build the 1.0 m line, use linear interpolation between the
c   model topo data.  Include the end point.
        do i=1,nx-1
          do ii=1,nnx
            iii = (i - 1) * nnx + ii
            x1 = 0.0
            x = real(ii - 1) * dx
            y2 = topo_line(i+1)
            y1 = topo_line(i)
            topo_1m(iii) = y1 + ((y2 - y1)/deltax) * (x - x1)
          enddo
        enddo
        topo_1m((nx - 1) * nnx + 1) = topo_line(nx)

c Use this topo array to be the starting point for generating the
c   Tabler surfaces.
        do i=1,maxlines
          tabler_1m(i) = topo_1m(i)
          drift_start_topo_1m(i) = 0.0
        enddo

c Run the Tabler model.
        do i=istart,iend,-1
          if (irotate_flag.eq.2) then
            t1 = tabler_1m(i)
            t2 = tabler_1m(i-31)
            t3 = tabler_1m(i-31-11)
            t4 = tabler_1m(i-31-21)
            t5 = tabler_1m(i-31-32)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(i-32) = max(topo_1m(i-32),
     &        tabler_1m(i-31) + y * slope_adjust * dx)

c Create an erosion area 'erosion_dist' upwind of the tabler surface.
            npts = nint(erosion_dist/1.41)
            if (tabler_1m(i-32).ne.topo_1m(i-32)) then
              do ii=1,npts
                if (tabler_1m(i-32+ii).eq.topo_1m(i-32+ii))
     &            drift_start_topo_1m(i-32+ii) = -8888.0
              enddo
            endif
          else
            t1 = tabler_1m(i)
            t2 = tabler_1m(i-44)
            t3 = tabler_1m(i-44-15)
            t4 = tabler_1m(i-44-30)
            t5 = tabler_1m(i-44-45)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(i-45) = max(topo_1m(i-45),
     &        tabler_1m(i-44) + y * slope_adjust * dx)

c Create an erosion area 'erosion_dist' upwind of the tabler surface.
            npts = nint(erosion_dist/1.0)
            if (tabler_1m(i-45).ne.topo_1m(i-45)) then
              do ii=1,npts
                if (tabler_1m(i-45+ii).eq.topo_1m(i-45+ii))
     &            drift_start_topo_1m(i-45+ii) = -8888.0
              enddo
            endif
          endif
        enddo

c Extract the profile at the model grid points.
        do i=1,nx
          ii = (i - 1) * nnx + 1
          tabler(i) = tabler_1m(ii)
          drift_start_topo(i) = drift_start_topo_1m(ii)
        enddo

c Use the 1-D arrays to fill in the 2-D tabler-surface array.
        do i=1,nx
          if (irotate_flag.eq.2) then
            JJ = j + i - nx
            if (JJ.ge.1 .and. JJ.le.ny) then
              tabler_ee(i,JJ) = tabler(i) + drift_start_topo(i)
            endif
          else
            tabler_ee(i,j) = tabler(i) + drift_start_topo(i)
          endif
        enddo

      enddo

c Convert snow_traps back to actual snow depths instead of
c   depth plus topography.
      do j=1,ny
        do i=1,nx
          tabler_ee(i,j) = tabler_ee(i,j) - topo_land(i,j)
        enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tabler_s(nx,ny,topo_land,tabler_ss,deltay,
     &  erosion_dist,irotate_flag,slope_adjust)

c This subroutine uses Tabler (1975) to define equilibrium profiles
c   for the topographic drift catchments, for the case of winds
c   from the SOUTH and SOUTHEAST.

c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
c   topographic catchments.  Proceedings of the 43rd Annual Western
c   Snow Conference, San Diego, California, 87-97.

      implicit none

      include 'snowmodel.inc'

      integer nx,ny,i,j,jstart,jend,II,irotate_flag,nny,nnx,
     &  jj,jjj,ny_max_tabler,maxlines
      integer npts
      real deltay,y,erosion_dist,slope_adjust,xmax_slope,dy,test,
     &  x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

      real topo_land(nx_max,ny_max)
      real tabler_ss(nx_max,ny_max)
      real tabler(ny_max)
      real topo_line(ny_max)
      real drift_start_topo(ny_max)

      parameter (ny_max_tabler=ny_max*1000)
      real topo_1m(ny_max_tabler)
      real tabler_1m(ny_max_tabler)
      real drift_start_topo_1m(ny_max_tabler)

c This program:
c   1) takes the coarse model topography and generates 1.0-m grid
c        increment topo lines;
c   2) uses those to generate the Tabler surface profiles; and
c   3) extracts the profiles at the coarse model grid cells.

c Fill the snow array with topo data.
      do j=1,ny
        do i=1,nx
          tabler_ss(i,j) = topo_land(i,j)
        enddo
      enddo

c Define the length of the i-looping, depending on whether there
c   is any rotation done to get 45 deg winds.
      if (irotate_flag.eq.2) then
        nnx = nx+ny-1
      else
        nnx = nx
      endif

c Required parameters.
      if (irotate_flag.eq.2) then
        dy = 1.41
        test = amod(deltay/1.41,dy/1.41)
      else
        dy = 1.0
        test = amod(deltay,dy)
      endif

      xmax_slope = -0.20

c This Tabler program has not been made general enough to deal
c   with deltax and deltay values that are not evenly divisible
c   by 1.0 (or 1.41 for diagonal profiles).
      if (abs(test).gt.1.0e10-5) then
        print *,'To generate the Tabler surfaces, deltax and deltay'
        print *,'  must be evenly divisible by 1.0 or 1.41.'
        print *,'  deltay = ',deltay
        stop
      endif

c Define the number of 1-m grid cells in each model grid cell, and
c   calculate how many of these are in the entire nx domain.
      nny = nint(deltay/dy)
      maxlines = (ny - 1) * nny + 1

c Define the starting and ending points of the line we are going to
c   work with.  This is done to make sure we are never looking
c   outside the data for numbers to work with.
      jstart = 1
      if (irotate_flag.eq.2) then
        jend = maxlines - (32+11+10+11)
      else
        jend = maxlines - (45+15+15+15)
      endif

c Extract the line we are going to work with.
      do i=1,nnx
        if (irotate_flag.eq.2) then
          do j=1,ny
            II = i - j + 1
            drift_start_topo(j) = 0.0
            if (II.le.0) then
              tabler(j) = tabler(j-1)
              topo_line(j) = tabler(j-1)
            elseif (II.gt.nx) then
              tabler(j) = tabler_ss(nx,i-nx+1)
              topo_line(j) = tabler_ss(nx,i-nx+1)
            else
              tabler(j) = tabler_ss(II,j)
              topo_line(j) = tabler_ss(II,j)
            endif
          enddo
        else
          do j=1,ny
            tabler(j) = tabler_ss(i,j)
            topo_line(j) = topo_land(i,j)
            drift_start_topo(j) = 0.0
          enddo
        endif

c To build the 1.0 m line, use linear interpolation between the
c   model topo data.  Include the end point.
        do j=1,ny-1
          do jj=1,nny
            jjj = (j - 1) * nny + jj
            x1 = 0.0
            x = real(jj - 1) * dy
            y2 = topo_line(j+1)
            y1 = topo_line(j)
            topo_1m(jjj) = y1 + ((y2 - y1)/deltay) * (x - x1)
          enddo
        enddo
        topo_1m((ny - 1) * nny + 1) = topo_line(ny)

c Use this topo array to be the starting point for generating the
c   Tabler surfaces.
        do j=1,maxlines
          tabler_1m(j) = topo_1m(j)
          drift_start_topo_1m(j) = 0.0
        enddo

c Run the Tabler model.
        do j=jstart,jend
          if (irotate_flag.eq.2) then
            t1 = tabler_1m(j)
            t2 = tabler_1m(j+31)
            t3 = tabler_1m(j+31+11)
            t4 = tabler_1m(j+31+21)
            t5 = tabler_1m(j+31+32)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(j+32) = max(topo_1m(j+32),
     &        tabler_1m(j+31) + y * slope_adjust * dy)

c Create an erosion area 'erosion_dist' upwind of the tabler surface.
            npts = nint(erosion_dist/1.41)
            if (tabler_1m(j+32).ne.topo_1m(j+32)) then
              do ii=1,npts
                if (tabler_1m(j+32-jj).eq.topo_1m(j+32-jj))
     &            drift_start_topo_1m(j+32-jj) = -8888.0
              enddo
            endif
          else
            t1 = tabler_1m(j)
            t2 = tabler_1m(j+44)
            t3 = tabler_1m(j+44+15)
            t4 = tabler_1m(j+44+30)
            t5 = tabler_1m(j+44+45)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(j+45) = max(topo_1m(j+45),
     &        tabler_1m(j+44) + y * slope_adjust * dy)

c Create an erosion area 'erosion_dist' upwind of the tabler surface.
            npts = nint(erosion_dist/1.0)
            if (tabler_1m(j+45).ne.topo_1m(j+45)) then
              do ii=1,npts
                if (tabler_1m(j+45-jj).eq.topo_1m(j+45-jj))
     &            drift_start_topo_1m(j+45-jj) = -8888.0
              enddo
            endif
          endif
        enddo

c Extract the profile at the model grid points.
        do j=1,ny
          jj = (j - 1) * nny + 1
          tabler(j) = tabler_1m(jj)
          drift_start_topo(j) = drift_start_topo_1m(jj)
        enddo

c Use the 1-D arrays to fill in the 2-D tabler-surface array.
        do j=1,ny
          if (irotate_flag.eq.2) then
            II = i - j + 1
            if (II.ge.1 .and. II.le.nx) then
              tabler_ss(II,j) = tabler(j) + drift_start_topo(j)
            endif
          else
            tabler_ss(i,j) = tabler(j) + drift_start_topo(j)
          endif
        enddo

      enddo

c Convert snow_traps back to actual snow depths instead of
c   depth plus topography.
      do j=1,ny
        do i=1,nx
          tabler_ss(i,j) = tabler_ss(i,j) - topo_land(i,j)
        enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tabler_n(nx,ny,topo_land,tabler_nn,deltay,
     &  erosion_dist,irotate_flag,slope_adjust)

c This subroutine uses Tabler (1975) to define equilibrium profiles
c   for the topographic drift catchments, for the case of winds
c   from the NORTH and NORTHWEST.

c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
c   topographic catchments.  Proceedings of the 43rd Annual Western
c   Snow Conference, San Diego, California, 87-97.

      implicit none

      include 'snowmodel.inc'

      integer nx,ny,i,j,jstart,jend,II,irotate_flag,nny,nnx,
     &  jj,jjj,ny_max_tabler,maxlines
      integer npts
      real deltay,y,erosion_dist,slope_adjust,xmax_slope,dy,test,
     &  x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

      real topo_land(nx_max,ny_max)
      real tabler_nn(nx_max,ny_max)
      real tabler(ny_max)
      real topo_line(ny_max)
      real drift_start_topo(ny_max)

      parameter (ny_max_tabler=ny_max*1000)
      real topo_1m(ny_max_tabler)
      real tabler_1m(ny_max_tabler)
      real drift_start_topo_1m(ny_max_tabler)

c This program:
c   1) takes the coarse model topography and generates 1.0-m grid
c        increment topo lines;
c   2) uses those to generate the Tabler surface profiles; and
c   3) extracts the profiles at the coarse model grid cells.

c Fill the snow array with topo data.
      do j=1,ny
        do i=1,nx
          tabler_nn(i,j) = topo_land(i,j)
        enddo
      enddo

c Define the length of the i-looping, depending on whether there
c   is any rotation done to get 45 deg winds.
      if (irotate_flag.eq.2) then
        nnx = nx+ny-1
      else
        nnx = nx
      endif

c Required parameters.
      if (irotate_flag.eq.2) then
        dy = 1.41
        test = amod(deltay/1.41,dy/1.41)
      else
        dy = 1.0
        test = amod(deltay,dy)
      endif

      xmax_slope = -0.20

c This Tabler program has not been made general enough to deal
c   with deltax and deltay values that are not evenly divisible
c   by 1.0 (or 1.41 for diagonal profiles).
      if (abs(test).gt.1.0e10-5) then
        print *,'To generate the Tabler surfaces, deltax and deltay'
        print *,'  must be evenly divisible by 1.0 or 1.41.'
        print *,'  deltay = ',deltay
        stop
      endif

c Define the number of 1-m grid cells in each model grid cell, and
c   calculate how many of these are in the entire nx domain.
      nny = nint(deltay/dy)
      maxlines = (ny - 1) * nny + 1

c Define the starting and ending points of the line we are going to
c   work with.  This is done to make sure we are never looking
c   outside the data for numbers to work with.
      jstart = maxlines
      if (irotate_flag.eq.2) then
        jend = 1 + (32+11+10+11)
      else
        jend = 1 + (45+15+15+15)
      endif

c Extract the line we are going to work with.
      do i=1,nnx
        if (irotate_flag.eq.2) then
          do j=1,ny
            II = i - j + 1
            drift_start_topo(j) = 0.0
            if (II.le.0) then
              tabler(j) = tabler(j-1)
              topo_line(j) = tabler(j-1)
            elseif (II.gt.nx) then
              tabler(j) = tabler_nn(nx,i-nx+1)
              topo_line(j) = tabler_nn(nx,i-nx+1)
            else
              tabler(j) = tabler_nn(II,j)
              topo_line(j) = tabler_nn(II,j)
            endif
          enddo
        else
          do j=1,ny
            tabler(j) = tabler_nn(i,j)
            topo_line(j) = topo_land(i,j)
            drift_start_topo(j) = 0.0
          enddo
        endif

c To build the 1.0 m line, use linear interpolation between the
c   model topo data.  Include the end point.
        do j=1,ny-1
          do jj=1,nny
            jjj = (j - 1) * nny + jj
            x1 = 0.0
            x = real(jj - 1) * dy
            y2 = topo_line(j+1)
            y1 = topo_line(j)
            topo_1m(jjj) = y1 + ((y2 - y1)/deltay) * (x - x1)
          enddo
        enddo
        topo_1m((ny - 1) * nny + 1) = topo_line(ny)

c Use this topo array to be the starting point for generating the
c   Tabler surfaces.
        do j=1,maxlines
          tabler_1m(j) = topo_1m(j)
          drift_start_topo_1m(j) = 0.0
        enddo

c Run the Tabler model.
        do j=jstart,jend,-1
          if (irotate_flag.eq.2) then
            t1 = tabler_1m(j)
            t2 = tabler_1m(j-31)
            t3 = tabler_1m(j-31-11)
            t4 = tabler_1m(j-31-21)
            t5 = tabler_1m(j-31-32)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(j-32) = max(topo_1m(j-32),
     &        tabler_1m(j-31) + y * slope_adjust * dy)

c Create an erosion area 'erosion_dist' upwind of the tabler surface.
            npts = nint(erosion_dist/1.41)
            if (tabler_1m(j-32).ne.topo_1m(j-32)) then
              do ii=1,npts
                if (tabler_1m(j-32+jj).eq.topo_1m(j-32+jj))
     &            drift_start_topo_1m(j-32+jj) = -8888.0
              enddo
            endif
          else
            t1 = tabler_1m(j)
            t2 = tabler_1m(j-44)
            t3 = tabler_1m(j-44-15)
            t4 = tabler_1m(j-44-30)
            t5 = tabler_1m(j-44-45)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(j-45) = max(topo_1m(j-45),
     &        tabler_1m(j-44) + y * slope_adjust * dy)

c Create an erosion area 'erosion_dist' upwind of the tabler surface.
            npts = nint(erosion_dist/1.0)
            if (tabler_1m(j-45).ne.topo_1m(j-45)) then
              do ii=1,npts
                if (tabler_1m(j-45+jj).eq.topo_1m(j-45+jj))
     &            drift_start_topo_1m(j-45+jj) = -8888.0
              enddo
            endif
          endif
        enddo

c Extract the profile at the model grid points.
        do j=1,ny
          jj = (j - 1) * nny + 1
          tabler(j) = tabler_1m(jj)
          drift_start_topo(j) = drift_start_topo_1m(jj)
        enddo

c Use the 1-D arrays to fill in the 2-D tabler-surface array.
        do j=1,ny
          if (irotate_flag.eq.2) then
            II = i - j + 1
            if (II.ge.1 .and. II.le.nx) then
              tabler_nn(II,j) = tabler(j) + drift_start_topo(j)
            endif
          else
            tabler_nn(i,j) = tabler(j) + drift_start_topo(j)
          endif
        enddo

      enddo

c Convert snow_traps back to actual snow depths instead of
c   depth plus topography.
      do j=1,ny
        do i=1,nx
          tabler_nn(i,j) = tabler_nn(i,j) - topo_land(i,j)
        enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine surface_snow_1(Tair,windspd,prec,ro_soft_snow,
     &    Utau_t,ro_soft_snow_old,dt,snow_z0,ht_windobs,
     &    ro_nsnow)

      implicit none

      real C,alfa,ro_min,ro_max,prec,Tair,ro_nsnow,dt,Tf,
     &  ro_offset,windspd,Utau_t,ro_soft_snow_old,ro_soft_snow,
     &  windspd_2m,snow_z0,ht_windobs

c Define the density rate coefficients.
      C = 0.10

c Define alfa. 
      alfa = 0.2

c Define the minimum and maximum snow density that should be
c   simulated.
      ro_min = 50.0
      ro_max = 450.0

c Freezing temperature.
      Tf = 273.16

c Calculate the new snow density.  First calculate this under the
c   assumption of no wind using the standard SnowModel formulation,
c   then calculate the offset for wind speeds > 5 m/s.
      if (prec.gt.0.0) then

c Calculate the 2-m wind speed.
        windspd_2m = windspd *
     &    log(2.0/snow_z0)/log(ht_windobs/snow_z0)

c To define the offset I assumed under cold conditions (ro_nsnow =
c   50.0 kg/m3):
c   1) 24-hour ave wind speed >= 20 m/s gives ro_nsnow = 350 kg/m3.
c   2) 24-hour ave wind speed = 5 m/s gives ro_nsnow = 150 kg/m3,
c      (actually what really matters here is the density difference
c      of 200 kg/m3 for the wind speed difference from 5 to 20 m/s).
c   3) 24-hour ave wind speed < 5 m/s gives ro_nsnow = ro_nsnow.
c   4) It is appropriate to use an exponential decay function
c      between the 5 and 20 m/s values.
        if (windspd_2m.lt.5.0) then
          ro_offset = 0.0
        else
          ro_offset = 25.0 +
     &      250.0 * (1.0 - exp(-(alfa*(windspd_2m - 5.0))))
        endif
        ro_soft_snow = ro_nsnow + ro_offset
        ro_soft_snow = min(ro_soft_snow,ro_max)
        ro_soft_snow = max(ro_soft_snow,ro_min)

        if (ro_soft_snow.le.300.0) then
          Utau_t = 0.10 * exp(0.003 * ro_soft_snow)
        else
          Utau_t = 0.005 * exp(0.013 * ro_soft_snow)
        endif
        ro_soft_snow_old = ro_soft_snow
      else
        call surface_snow_2(ro_soft_snow_old,ro_soft_snow,Utau_t,
     &    dt,Tair,windspd_2m,C,ro_max,ro_min,alfa,Tf)
        ro_soft_snow_old = ro_soft_snow
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine surface_snow_2(ro_soft_snow_old,ro_soft_snow,Utau_t,
     &  dt,Tair,windspd_2m,C,ro_max,ro_min,alfa,Tf)

      implicit none

      real Tf,Tsnow,Tair,ro_soft_snow,dt,A1,A2,B,U,C,
     &  windspd_2m,Utau_t,ro_soft_snow_old,ro_max,ro_min,alfa

      A1 = 0.0013
      A2 = 0.021
      B = 0.08

c Evolve the near-surface snow density under the influence of
c   temperature and snow-transporting wind speeds.

c Assume that the near-surface snow temperature equals the air
c   temperature, but is not above the melting point.
      Tsnow = min(Tf,Tair)

c Update the snow density of the soft snow layer.  Eliminate the
c   wind speed influence for speeds below 5 m/s, but account for it
c   if speeds are >= 5 m/s.
      if (windspd_2m.ge.5.0) then
        U = 5.0 + 15.0 * (1.0 - exp(-(alfa*(windspd_2m - 5.0))))
      else
        U = 1.0
      endif

      ro_soft_snow = ro_soft_snow_old + dt *
     &  (C * A1 * U * ro_soft_snow_old *
     &  exp((- B)*(Tf-Tsnow)) * exp((- A2)*ro_soft_snow_old))

c Bound the calculated density.
      ro_soft_snow = min(ro_max,ro_soft_snow)
      ro_soft_snow = max(ro_min,ro_soft_snow)

c Calculate the snow threshold friction velocity.
      if (ro_soft_snow.le.300.0) then
        Utau_t = 0.10 * exp(0.003 * ro_soft_snow)
      else
        Utau_t = 0.005 * exp(0.013 * ro_soft_snow)
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

