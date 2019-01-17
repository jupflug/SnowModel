ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c snowmodel_main.f
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c All units are in m, kg, s, K.
c
c The author of this code is:
c
c   Dr. Glen E. Liston
c   Cooperative Institute for Research
c     in the Atmosphere (CIRA)         |
c   Colorado State University
c   Fort Collins, Colorado 80523-1375
c
c   Voice: (970) 491-8220
c   FAX: (970) 491-8241
c
c   glen.liston@colostate.edu

      implicit none

      include 'snowmodel.inc'

      include 'snowmodel_vars.inc'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc  INITIALIZE THE MODEL  cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Read the input parameters.
      CALL READPARAM_CODE(dt,deltax,deltay,Utau_t_flag,
     &  subgrid_flag,erosion_dist,tp_scale,twolayer_flag,
     &  bc_flag,curve_len_scale,slopewt,curvewt,ht_windobs,
     &  ht_rhobs,ro_snow,snow_d_init_const,const_veg_flag,
     &  vegsnowdepth,nx,ny,max_iter,met_input_fname,xmn,ymn,
     &  iyear_init,imonth_init,iday_init,xhour_init,undef,ifill,
     &  iobsint,dn,xlat,i_tair_flag,i_rh_flag,i_wind_flag,
     &  i_solar_flag,i_prec_flag,isingle_stn_flag,igrads_metfile,
     &  windspd_min,icond_flag,run_micromet,run_enbal,run_snowpack,
     &  run_snowtran,topoflag,topoveg_fname,snowtran_output_fname,
     &  micromet_output_fname,enbal_output_fname,Utau_t_const,
     &  snowpack_output_fname,print_micromet,print_enbal,
     &  print_snowpack,print_snowtran,i_longwave_flag,print_user,
     &  iprint_inc,ascii_topoveg,topo_ascii_fname,veg_ascii_fname,
     &  irun_data_assim,lapse_rate_user_flag,
     &  iprecip_lapse_rate_user_flag,use_shortwave_obs,
     &  use_longwave_obs,use_sfc_pressure_obs,calc_subcanopy_met,
     &  sfc_sublim_flag,gap_frac,cloud_frac_factor,
     &  albedo_glacier,barnes_lg_domain,n_stns_used,tabler_dir,
     &  slope_adjust,lat_solar_flag,UTC_flag,iveg_ht_flag,
     &  ihrestart_flag,ihrestart_inc,i_dataassim_loop,tsls_threshold,
     &  dz_snow_min,print_multilayer,multilayer_snowpack,max_layers,
     &  multilayer_output_fname,izero_snow_date,curve_lg_scale_flag,
     &  check_met_data,seaice_run,snowmodel_line_flag,wind_lapse_rate,
     &  albedo_diff,al_max,al_min,al_dec_cold,al_dec_melt,
     &  fc_param,t_avg)

c This loop runs the correction/data assimilation adjustment
c   iterations.
      if (ihrestart_flag.ge.0) then
        if (i_dataassim_loop.lt.0.0) then
          i_corr_start = 2
        else
          i_corr_start = 1
        endif
      else
        i_corr_start = 1
      endif

      do icorr_factor_loop=i_corr_start,irun_data_assim+1

c Perform the correction (precipitation and melt) factor
c   calculations.
        if (irun_data_assim.eq.1 .and. icorr_factor_loop.eq.2) then
          CALL DATAASSIM_USER(nx,ny,icorr_factor_index,
     &      corr_factor,max_iter,deltax,deltay,xmn,ymn,nobs_dates,
     &      print_inc,iday_init,imonth_init,iyear_init)
          if (ihrestart_flag.ge.-1) then
            CALL HRESTART_SAVE_DA(nx,ny,max_iter,corr_factor,
     &        icorr_factor_index,nobs_dates)
          endif
        endif

c Perform a variety of preprocessing and model setup steps, like
c   read in topography and vegetation arrays, open input and output
c   files, etc.
        CALL PREPROCESS_CODE(topoveg_fname,const_veg_flag,
     &    vegtype,veg_z0,vegsnowdepth,fetch,xmu,C_z,h_const,
     &    wind_min,Up_const,dz_susp,ztop_susp,fall_vel,Ur_const,
     &    ro_water,ro_air,gravity,vonKarman,pi,twopio360,snow_z0,
     &    nx,ny,sum_sprec,sum_qsubl,sum_trans,sum_unload,topo,
     &    topo_land,snow_d,topoflag,snow_d_init,snow_d_init_const,
     &    soft_snow_d,met_input_fname,igrads_metfile,deltax,deltay,
     &    snowtran_output_fname,micromet_output_fname,
     &    enbal_output_fname,snowpack_output_fname,print_micromet,
     &    print_enbal,print_snowpack,print_snowtran,run_micromet,
     &    run_enbal,run_snowpack,run_snowtran,ro_snow_grid,swe_depth,
     &    sum_runoff,sum_prec,ro_snow,twolayer_flag,sum_Qcs,
     &    canopy_int,ascii_topoveg,topo_ascii_fname,icorr_factor_loop,
     &    veg_ascii_fname,undef,isingle_stn_flag,max_iter,
     &    i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag,sum_glacmelt,
     &    snow_depth,sum_d_canopy_int,corr_factor,icorr_factor_index,
     &    sum_sfcsublim,barnes_lg_domain,n_stns_used,k_stn,xmn,ymn,
     &    ro_soft_snow_old,sum_swemelt,xlat,lat_solar_flag,xlat_grid,
     &    xlon_grid,UTC_flag,dt,swe_depth_old,canopy_int_old,
     &    vegsnowd_xy,iveg_ht_flag,ihrestart_flag,i_dataassim_loop,
     &    multilayer_snowpack,max_layers,multilayer_output_fname,
     &    print_multilayer,JJ,tslsnowfall,tsls_threshold,
     &    irun_data_assim,izero_snow_date,iclear_mn,iclear_dy,
     &    xclear_hr,dy_snow,swe_lyr,ro_layer,T_old,gamma,icond_flag,
     &    curve_lg_scale_flag,curve_wt_lg,check_met_data,seaice_run,
     &    snowmodel_line_flag,xg_line,yg_line,print_user)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccc  RUN THE MODEL  cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Loop through the number of time steps (iterations) in the
c   simulation.
        if (ihrestart_flag.ge.0) then
          iter_start = ihrestart_flag + 1
        else
          iter_start = 1
        endif

        do iter=iter_start,max_iter

c Distribute the meteorological station data.
          if (run_micromet.eq.1.0) then
            CALL MICROMET_CODE(nx,ny,xmn,ymn,deltax,deltay,
     &        iyear_init,imonth_init,iday_init,xhour_init,dt,undef,
     &        ifill,iobsint,dn,iter,curve_len_scale,slopewt,curvewt,
     &        topo,curvature,terrain_slope,slope_az,Tair_grid,
     &        rh_grid,uwind_grid,vwind_grid,Qsi_grid,prec_grid,
     &        i_tair_flag,i_rh_flag,i_wind_flag,i_solar_flag,
     &        i_prec_flag,isingle_stn_flag,igrads_metfile,
     &        windspd_grid,winddir_grid,windspd_flag,winddir_flag,
     &        sprec,windspd_min,Qli_grid,i_longwave_flag,vegtype,
     &        forest_LAI,iyear,imonth,iday,xhour,corr_factor,
     &        icorr_factor_index,lapse_rate_user_flag,
     &        iprecip_lapse_rate_user_flag,use_shortwave_obs,
     &        use_longwave_obs,use_sfc_pressure_obs,sfc_pressure,
     &        run_enbal,run_snowpack,calc_subcanopy_met,vegsnowd_xy,
     &        gap_frac,cloud_frac_factor,barnes_lg_domain,n_stns_used,
     &        k_stn,xlat_grid,xlon_grid,UTC_flag,icorr_factor_loop,
     &        snowmodel_line_flag,xg_line,yg_line,irun_data_assim,
     &        wind_lapse_rate,prec_grid_sol)
            
            if (print_micromet.eq.1.0) then
              if (mod(iter,iprint_inc).eq.0) then
                write(81,rec=iter/iprint_inc)
     &            ((Tair_grid(i,j)-273.16,i=1,nx),j=1,ny),
     &            ((rh_grid(i,j),i=1,nx),j=1,ny),
     &            ((windspd_grid(i,j),i=1,nx),j=1,ny),
     &            ((winddir_grid(i,j),i=1,nx),j=1,ny),
     &            ((Qsi_grid(i,j),i=1,nx),j=1,ny),
     &            ((Qli_grid(i,j),i=1,nx),j=1,ny),
     &            ((prec_grid(i,j),i=1,nx),j=1,ny),
     &            ((sprec(i,j),i=1,nx),j=1,ny)
              endif
            endif
          endif

c Perform a surface energy balance over the domain.
          if (run_enbal.eq.1.0) then
            CALL ENBAL_CODE(nx,ny,Tair_grid,uwind_grid,sfc_pressure,
     &        vwind_grid,rh_grid,Tsfc,Qsi_grid,Qli_grid,Qle,Qh,Qe,
     &        Qc,Qm,e_balance,Qf,snow_d,ht_windobs,icond_flag,
     &        albedo,snow_z0,veg_z0,vegtype,undef,albedo_glacier,
     &        dy_snow,T_old,gamma,JJ,prec_grid,albedo_diff,al_max,
     &        al_min,al_dec_cold,al_dec_melt,dt,sprec,Qp)

            if (print_enbal.eq.1.0) then
              if (mod(iter,iprint_inc).eq.0) then
                write(82,rec=iter/iprint_inc)
     &            ((Tsfc(i,j)-273.16,i=1,nx),j=1,ny),
     &            ((Qsi_grid(i,j),i=1,nx),j=1,ny),
     &            ((Qli_grid(i,j),i=1,nx),j=1,ny),
     &            ((Qle(i,j),i=1,nx),j=1,ny),
     &            ((Qh(i,j),i=1,nx),j=1,ny),
     &            ((Qe(i,j),i=1,nx),j=1,ny),
     &            ((Qc(i,j),i=1,nx),j=1,ny),
     &            ((Qm(i,j),i=1,nx),j=1,ny),
     &            ((albedo(i,j),i=1,nx),j=1,ny),
     &            ((e_balance(i,j),i=1,nx),j=1,ny)
              endif
            endif
          endif

c Evolve the snowpack according to the defined melt and
c   precipitation inputs.
          if (run_snowpack.eq.1.0) then
            CALL SNOWPACK_CODE(nx,ny,Tair_grid,rh_grid,ro_nsnow,
     &        dt,swe_depth,Tsfc,snow_d,prec_grid,runoff,Qm,rain,
     &        sprec,iter,w_balance,sum_prec,sum_runoff,xro_snow,
     &        undef,ro_snow,ro_snow_grid,soft_snow_d,sum_sprec,
     &        snow_depth,windspd_grid,Qsi_grid,sum_Qcs,canopy_int,
     &        Qcs,vegtype,forest_LAI,albedo,glacier_melt,
     &        canopy_unload,sum_unload,sum_glacmelt,run_snowtran,
     &        swemelt,d_canopy_int,sum_d_canopy_int,snow_d_init,
     &        sfc_pressure,Qe,sfc_sublim_flag,sum_sfcsublim,
     &        sum_swemelt,corr_factor,icorr_factor_index,swesublim,
     &        swe_depth_old,canopy_int_old,JJ,max_layers,melt_flag,
     &        ro_snowmax,tsls_threshold,dz_snow_min,tslsnowfall,
     &        change_layer,dy_snow,swe_lyr,ro_layer,T_old,gamma,
     &        multilayer_snowpack,seaice_run,seaice_conc,
     &        fc_param,t_avg,Saturn)
          endif

c Run the blowing-snow model.
          if (run_snowtran.eq.1.0) then
            CALL SNOWTRAN_CODE(bc_flag,bs_flag,C_z,
     &        conc_salt,deltax,deltay,dh_salt,dh_salt_u,dh_salt_v,
     &        dh_susp,dh_susp_u,dh_susp_v,dt,dz_susp,fall_vel,fetch,
     &        gravity,h_const,h_star,ht_rhobs,ht_windobs,index_ue,
     &        index_uw,index_vn,index_vs,iter,nx,ny,pi,Qsalt,Qsalt_max,
     &        Qsalt_maxu,Qsalt_maxv,Qsalt_u,Qsalt_v,Qsubl,Qsusp,
     &        Qsusp_u,Qsusp_v,rh_grid,ro_air,ro_snow,ro_water,snow_d,
     &        snow_d_init,snow_z0,soft_snow_d,sprec,sum_glacmelt,
     &        subgrid_flag,wbal_salt,wbal_susp,wbal_qsubl,sum_sprec,
     &        tabler_ee,tabler_ne,tabler_nn,tabler_nw,tabler_se,
     &        tabler_ss,tabler_sw,tabler_ww,tair_grid,topo,topo_land,
     &        topoflag,tp_scale,twolayer_flag,Up_const,Ur_const,Utau,
     &        Utau_t,uwind_grid,veg_z0,vegsnowd_xy,vegtype,vonKarman,
     &        vwind_grid,wind_min,winddir_flag,winddir_grid,
     &        windspd_flag,windspd_grid,xmu,z_0,ztop_susp,erosion_dist,
     &        run_enbal,run_snowpack,wbal_subgrid,sum_qsubl,sum_trans,
     &        swe_depth,snow_depth,ro_snow_grid,sum_prec,sum_runoff,
     &        sum_Qcs,canopy_int,w_balance,sum_sfcsublim,tabler_dir,
     &        slope_adjust,Utau_t_const,Utau_t_flag,ro_soft_snow_old,
     &        ro_soft_snow,ro_nsnow,prec_grid,Qcs,runoff,d_canopy_int,
     &        glacier_melt,swe_depth_old,swesublim,canopy_unload,
     &        canopy_int_old,iter_start,multilayer_snowpack,swe_lyr,
     &        JJ,dy_snow,ro_layer,curve_lg_scale_flag,curve_wt_lg,
     &        seaice_run,seaice_conc,tslsnowfall,T_old,tsls_threshold)
          endif

c Save the outputs from the SNOWPACK and SNOWTRAN routines.
          if (run_snowpack.eq.1.0 .and. print_snowpack.eq.1.0) then
            if (mod(iter,iprint_inc).eq.0) then
              write(83,rec=iter/iprint_inc)
     &          ((snow_depth(i,j),i=1,nx),j=1,ny),
     &          ((xro_snow(i,j),i=1,nx),j=1,ny),
     &          ((swe_depth(i,j),i=1,nx),j=1,ny),
     &          ((runoff(i,j),i=1,nx),j=1,ny),
     &          ((rain(i,j),i=1,nx),j=1,ny),
     &          ((sprec(i,j),i=1,nx),j=1,ny),
     &          ((Qcs(i,j),i=1,nx),j=1,ny),
     &          ((Tsfc(i,j),i=1,nx),j=1,ny),
     &          ((sum_Qcs(i,j),i=1,nx),j=1,ny),
     &          ((sum_prec(i,j),i=1,nx),j=1,ny),
     &          ((sum_sprec(i,j),i=1,nx),j=1,ny),
     &          ((sum_unload(i,j),i=1,nx),j=1,ny),
     &          ((sum_runoff(i,j),i=1,nx),j=1,ny),
     &          ((sum_swemelt(i,j),i=1,nx),j=1,ny),
     &          ((sum_sfcsublim(i,j),i=1,nx),j=1,ny),
     &          ((w_balance(i,j),i=1,nx),j=1,ny)
            endif
          endif

          if (run_snowtran.eq.1.0 .and. print_snowtran.eq.1.0) then
            if (mod(iter,iprint_inc).eq.0) then
              write(84,rec=iter/iprint_inc)
     &          ((snow_d(i,j),i=1,nx),j=1,ny),
     &          ((wbal_qsubl(i,j),i=1,nx),j=1,ny),
     &          ((wbal_salt(i,j),i=1,nx),j=1,ny),
     &          ((wbal_susp(i,j),i=1,nx),j=1,ny),
     &          ((wbal_subgrid(i,j),i=1,nx),j=1,ny),
     &          ((sum_qsubl(i,j),i=1,nx),j=1,ny),
     &          ((sum_trans(i,j),i=1,nx),j=1,ny)
            endif
          endif

c Note that here we write out the entire potential vertical domain
c   (nz_max), because it may have been filled with snow at some point
c   during the simulation.
          if (run_snowpack.eq.1.0 .and. multilayer_snowpack.eq.1 .and.
     &      print_multilayer.eq.1.0) then
            if (mod(iter,iprint_inc).eq.0) then
              write(85,rec=iter/iprint_inc)
     &          ((snow_depth(i,j),i=1,nx),j=1,ny),
     &          ((swe_depth(i,j),i=1,nx),j=1,ny),
     &          ((real(JJ(i,j)),i=1,nx),j=1,ny),
     &          ((change_layer(i,j),i=1,nx),j=1,ny),
     &          (((T_old(i,j,k),i=1,nx),j=1,ny),k=1,nz_max),
c     &          (((gamma(i,j,k),i=1,nx),j=1,ny),k=1,nz_max),
     &          (((Saturn(i,j,k),i=1,nx),j=1,ny),k=1,nz_max),
     &          (((ro_layer(i,j,k),i=1,nx),j=1,ny),k=1,nz_max),
     &          (((swe_lyr(i,j,k),i=1,nx),j=1,ny),k=1,nz_max),
     &          (((dy_snow(i,j,k),i=1,nx),j=1,ny),k=1,nz_max),
     &          (((real(melt_flag(i,j,k)),i=1,nx),j=1,ny),k=1,nz_max)
            endif
          endif

c The call to outputs_user is available to provide user-defined
c   outputs.  These might be special-case situations, like just
c   writing out data at the end of every day, writing out a few
c   grid cells, saving each data arrays to individual files, etc.
          if (print_user.eq.1.0) then
            CALL OUTPUTS_USER(nx,ny,iter,Tair_grid,rh_grid,
     &        uwind_grid,vwind_grid,windspd_grid,winddir_grid,
     &        Qsi_grid,Qli_grid,prec_grid,Tsfc,Qle,Qh,Qe,Qc,Qm,Qf,
     &        e_balance,snow_depth,xro_snow,swe_depth,ro_nsnow,
     &        runoff,rain,sprec,sum_prec,sum_runoff,w_balance,
     &        snow_d,topo_land,wbal_qsubl,sum_sprec,wbal_salt,
     &        wbal_susp,ro_snow_grid,sum_Qcs,canopy_int,Qcs,
     &        iyear,imonth,iday,xhour,undef,deltax,xmn,ymn,
     &        wbal_subgrid,canopy_unload,sum_qsubl,sum_trans,
     &        sum_unload,sum_glacmelt,glacier_melt,swemelt,
     &        iprint_inc,sfc_pressure,sum_swemelt,albedo,
     &        icorr_factor_loop,swesublim,vegtype,iter_start,
     &        seaice_run,print_inc)
          endif

c For multi-year simulations, sometimes it is desirable to zero
c   out the snow cover arrays on a certain summer date, to prevent
c   glaciers from forming.
          if (imonth.eq.iclear_mn .and. iday.eq.iclear_dy .and.
     &      xhour.eq.xclear_hr) then
            CALL ZERO_SNOW(nx,ny,snow_depth,ro_snow_grid,ro_snow,
     &        swe_depth,swe_depth_old,canopy_int_old,JJ,sum_swemelt,
     &        tslsnowfall,dy_snow,swe_lyr,ro_layer,T_old,sum_sprec,
     &        multilayer_snowpack,tsls_threshold)
          endif

c Save the history restart information.
          if (ihrestart_flag.ge.-1) then
            if (mod(iter,ihrestart_inc).eq.0
     &        .or. iter.eq.max_iter) then
              CALL HRESTART_SAVE(nx,ny,iter,snow_d,snow_depth,
     &          canopy_int,soft_snow_d,ro_snow_grid,swe_depth,
     &          ro_soft_snow_old,snow_d_init,swe_depth_old,
     &          canopy_int_old,topo,sum_sprec,icorr_factor_loop,
     &          max_iter)
            endif
          endif

        enddo

      enddo

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

