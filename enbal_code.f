c enbal_code.f

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ENBAL_CODE(nx,ny,Tair_grid,uwind_grid,sfc_pressure,
     &  vwind_grid,rh_grid,Tsfc,Qsi_grid,Qli_grid,Qle,Qh,Qe,
     &  Qc,Qm,e_balance,Qf,snow_d,ht_windobs,icond_flag,
     &  albedo,snow_z0,veg_z0,vegtype,undef,albedo_glacier,dy_snow,
     &  T_old,gamma,JJ,prec_grid,albedo_diff,al_max,al_min,al_dec_cold,
     &  al_dec_melt,dt,sprec,Qp,albedo_flag)

      implicit none

      include 'snowmodel.inc'

      real Tair_grid(nx_max,ny_max)
      real rh_grid(nx_max,ny_max)
      real uwind_grid(nx_max,ny_max)
      real vwind_grid(nx_max,ny_max)
      real Qsi_grid(nx_max,ny_max)
      real Qli_grid(nx_max,ny_max)
      real albedo(nx_max,ny_max)
      real vegtype(nx_max,ny_max)
      real veg_z0(nx_max,ny_max)
      real prec_grid(nx_max,ny_max)
      real sprec(nx_max,ny_max)
      real Qp(nx_max,ny_max)

      real Tsfc(nx_max,ny_max),Qle(nx_max,ny_max),
     &  Qh(nx_max,ny_max),Qe(nx_max,ny_max),Qc(nx_max,ny_max),
     &  Qm(nx_max,ny_max),e_balance(nx_max,ny_max),Qf(nx_max,ny_max),
     &  snow_d(nx_max,ny_max),sfc_pressure(nx_max,ny_max)

      real snow_z0,veg_z0_tmp,windspd,ht_windobs,undef,
     &  albedo_glacier,count_Tsfc_not_converged,albedo_diff,al_max,
     &  al_min,al_dec_cold,al_dec_melt,dt

      integer i,j,nx,ny,icond_flag,k,albedo_flag

      integer JJ(nx_max,ny_max)
      real dy_snow(nx_max,ny_max,nz_max)
      real T_old(nx_max,ny_max,nz_max)
      real gamma(nx_max,ny_max,nz_max)
      real dy_snow_z(2)
      real T_old_z(2)
      real gamma_z(2)

      print *,'   solving the energy balance'
      count_Tsfc_not_converged = 0.0

      do j=1,ny
        do i=1,nx

          windspd = sqrt(uwind_grid(i,j)**2+vwind_grid(i,j)**2)
          veg_z0_tmp = veg_z0(i,j)

c Extract the vertical column for this i,j point, and send it
c   to the subroutine. *** Note that I should use f95, then I would
c   not have to do this (I could pass in subsections of the arrays).
          if (icond_flag.eq.1) then
            do k=1,2
              dy_snow_z(k) = dy_snow(i,j,k)
              T_old_z(k) = T_old(i,j,k)
              gamma_z(k) = gamma(i,j,k)
            enddo
          endif

          CALL ENBAL_CORE(Tair_grid(i,j),windspd,rh_grid(i,j),
     &      Tsfc(i,j),Qsi_grid(i,j),Qli_grid(i,j),Qle(i,j),Qh(i,j),
     &      Qe(i,j),Qc(i,j),Qm(i,j),e_balance(i,j),Qf(i,j),undef,
     &      sfc_pressure(i,j),snow_d(i,j),ht_windobs,
     &      icond_flag,albedo(i,j),snow_z0,veg_z0_tmp,vegtype(i,j),
     &      albedo_glacier,dy_snow_z,T_old_z,gamma_z,JJ(i,j),
     &      count_Tsfc_not_converged,prec_grid(i,j),albedo_diff,al_max,
     &      al_min,al_dec_cold,al_dec_melt,dt,sprec(i,j),Qp(i,j),
     &      albedo_flag)
        enddo
      enddo

c Calculate the % of the grid cells that did not converge during
c   this time step.
      count_Tsfc_not_converged = 100.0* count_Tsfc_not_converged /
     &  real(nx*ny)

c Set the not-converged threshold to be 1% of the grid cells.
      if (count_Tsfc_not_converged.gt.1.0) then
        print *,'Over 1% of the grid cells failed to converge'
        print *,'  in the Tsfc energy balance calculation. This'
        print *,'  usually means there in a problem with the'
        print *,'  atmopheric forcing inputs, or that windspd_min'
        print *,'  in snowmodel.par is set too low; like less than'
        print *,'  1 m/s.'
        print *
        print *,'% Tsfc not converged = ',count_Tsfc_not_converged

      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ENBAL_CORE(Tair,windspd,rh,
     &    Tsfc,Qsi,Qli,Qle,Qh,
     &    Qe,Qc,Qm,e_balance,Qf,undef,
     &    sfc_pressure,snow_d,ht_windobs,
     &    icond_flag,albedo,snow_z0,veg_z0_tmp,vegtype,
     &    albedo_glacier,dy_snow_z,T_old_z,gamma_z,JJ,
     &    count_Tsfc_not_converged,prec,albedo_diff,al_max,al_min,
     &    al_dec_cold,al_dec_melt,dt,sprec,Qp,albedo_flag)

c This is the FORTRAN code which implements the surface energy
c   balance model used in:
c
c   (1) "Local Advection of Momentum, Heat, and Moisture during the
c       Melt of Patchy Snow Covers", G. E. Liston, J. Applied 
c       Meteorology, 34, 1705-1715, 1995.
c
c   (2) "An Energy-Balance Model of Lake-Ice Evolution", G. E.
c       Liston, D. K. Hall, J. of Glaciology, (41), 373-382, 1995.
c
c   (3) "Sensitivity of Lake Freeze-Up and Break-Up to Climate
c       Change: A Physically Based Modeling Study", G. E. Liston,
c       D. K. Hall, Annals of Glaciology, (21), 387-393, 1995.
c
c   (4) "Below-Surface Ice Melt on the Coastal Antarctic Ice
c       Sheet", G. E. Liston, and 4 others, J. of Glaciology, (45),
c       273-285, 1999.
c
c This version runs at sub-hourly to daily time steps.
c
c In addition, this version includes the influence of direct and
c   diffuse solar radiation, and the influence of topographic 
c   slope and aspect on incoming solar radiation.

      implicit none

      real Tair,windspd,rh,Tsfc,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,e_balance,
     &  Qf,sfc_pressure,snow_d,ht_windobs,albedo,snow_z0,
     &  veg_z0_tmp,vegtype,emiss_sfc,Stef_Boltz,ro_air,Cp,gravity,
     &  xls,xkappa,xLf,Tf,ro_water,Cp_water,ro_ice,z_0,ea,de_h,
     &  stability,es0,undef,albedo_glacier,count_Tsfc_not_converged,
     &  prec,albedo_diff,al_max,al_min,al_dec_cold,al_dec_melt,dt,
     &  sprec,Cp_snow,Qp,Twb

      integer icond_flag,JJ,albedo_flag

      real dy_snow_z(2)
      real T_old_z(2)
      real gamma_z(2)

c Define the constants used in the computations.
        CALL CONSTS_ENBAL(Stef_Boltz,ro_air,Cp,gravity,
     &    xls,xkappa,xLf,Tf,ro_water,Cp_water,ro_ice,emiss_sfc,
     &    Cp_snow)

c Define the surface characteristics based on the snow conditions.
        CALL GET_SFC(snow_d,albedo,z_0,Tf,Tair,snow_z0,veg_z0_tmp,
     &    vegtype,albedo_glacier,prec,albedo_diff,al_max,al_min,
     &    al_dec_cold,al_dec_melt,dt,sprec,albedo_flag)

c Atmospheric vapor pressure from relative humidity data.
        CALL VAPPRESS(ea,rh,Tair,Tf)

c Compute the turbulent exchange coefficients.
        CALL EXCOEFS(De_h,z_0,ht_windobs,windspd,xkappa)

c Compute the flux contribution due to conduction.
        CALL CONDUCT(icond_flag,Qc,dy_snow_z,T_old_z,gamma_z,JJ)

c J.PFLUG
c solve for the temperature of the precipitation (assuming prectemp = 
c the wet-bulb temperature) for advection.
        CALL SOLVEWB_EB(Twb,Tf,Tair,rh,xLs,Cp,sfc_pressure)
c END J.PFLUG

c Solve the energy balance for the surface temperature.
        CALL SFCTEMP(Tsfc,Tair,Qsi,Qli,ea,albedo,De_h,
     &    sfc_pressure,ht_windobs,windspd,ro_air,Cp,emiss_sfc,
     &    Stef_Boltz,gravity,xLs,xkappa,z_0,Tf,Qc,
     &    count_Tsfc_not_converged,prec,sprec,Cp_snow,Cp_water,xLf,
     &    Qp,dt,Twb)

c Make sure the snow surface temperature is <= 0 C.
        CALL MELTTEMP(Tsfc,Tf,snow_d,vegtype)

c Compute the stability function.
        CALL STABLEFN(stability,Tair,Tsfc,windspd,ht_windobs,
     &    gravity,xkappa,z_0)

c Compute the water vapor pressure at the surface.
        CALL VAPOR(es0,Tsfc,Tf)

c Compute the latent heat flux.
        CALL LATENT(Qe,De_h,stability,ea,es0,ro_air,xLs,
     &    sfc_pressure,undef,snow_d,vegtype)

c Compute the sensible heat flux.
        CALL SENSIBLE(Qh,De_h,stability,Tair,Tsfc,ro_air,Cp,
     &    undef,snow_d,vegtype)

c Compute the longwave flux emitted by the surface.
        CALL LONGOUT(Qle,Tsfc,emiss_sfc,Stef_Boltz)

c Compute the energy flux available for melting or freezing.
        CALL MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,Qf,Tsfc,
     &    Tf,Tair,windspd,ht_windobs,gravity,De_h,ea,ro_air,xLs,
     &    sfc_pressure,Cp,emiss_sfc,Stef_Boltz,snow_d,xkappa,z_0,
     &    vegtype,icond_flag,undef,dy_snow_z,T_old_z,gamma_z,JJ,Qp)

c Perform an energy balance check.
        CALL ENBAL(e_balance,albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,
     &    undef,snow_d,vegtype,Qp)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE VAPPRESS(ea,rh,Tair,Tf)

      implicit none

      real A,B,C,ea,rh,Tair,Tf

c Also see the VAPOR subroutine.

c Coeffs for saturation vapor pressure over water (Buck 1981).
c   Note: temperatures for Buck`s equations are in deg C, and
c   vapor pressures are in mb.  Do the adjustments so that the
c   calculations are done with temperatures in K, and vapor
c   pressures in Pa.

c Because I am interested in sublimation over snow during the
c   winter, do these calculations over ice.

c Over water.
c       A = 6.1121 * 100.0
c       B = 17.502
c       C = 240.97
c Over ice.
        A = 6.1115 * 100.0
        B = 22.452
        C = 272.55

c Atmospheric vapor pressure from relative humidity data.
      ea = rh / 100.0 * A * exp((B * (Tair - Tf))/(C + (Tair - Tf)))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,Qf,Tsfc,
     &  Tf,Tair,windspd,ht_windobs,gravity,De_h,ea,ro_air,xLs,
     &  sfc_pressure,Cp,emiss_sfc,Stef_Boltz,snow_d,xkappa,z_0,
     &  vegtype,icond_flag,undef,dy_snow_z,T_old_z,gamma_z,JJ,Qp)

      implicit none

      real albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,Qf,Tsfc,Tf,Tair,windspd,
     &  ht_windobs,gravity,De_h,ea,ro_air,xLs,sfc_pressure,Cp,
     &  emiss_sfc,Stef_Boltz,snow_d,xkappa,z_0,vegtype,xTsfc,
     &  xstability,xes0,xQe,xqh,xQle,xQc,undef,Qp

      integer icond_flag,JJ

      real dy_snow_z(2)
      real T_old_z(2)
      real gamma_z(2)

c If Qm is > 0, then this is the energy available for melting.
c   If Qm is < 0, then this is the energy available for freezing
c   liquid water in the snowpack.
c J.PFLUG
c heat advected by precipitation is also included
      if (snow_d.gt.0.0 .and. Tsfc.eq.Tf) then
c        Qm = (1.0-albedo) * Qsi + Qli + Qle + Qh + Qe + Qc + Qp
        Qm = (1.0-albedo) * Qsi + Qli + Qle + Qh + Qe + Qc
      elseif (vegtype.eq.20.0 .and. Tsfc.eq.Tf) then
c        Qm = (1.0-albedo) * Qsi + Qli + Qle + Qh + Qe + Qc + Qp
        Qm = (1.0-albedo) * Qsi + Qli + Qle + Qh + Qe + Qc
      else
        Qm = 0.0
      endif

      if (Tsfc.lt.Tf) then
        xTsfc = Tf
        CALL STABLEFN(xstability,Tair,xTsfc,windspd,ht_windobs,
     &    gravity,xkappa,z_0)
        CALL VAPOR(xes0,xTsfc,Tf)
        CALL LATENT(xQe,De_h,xstability,ea,xes0,ro_air,xLs,
     &    sfc_pressure,undef,snow_d,vegtype)
        CALL SENSIBLE(xQh,De_h,xstability,Tair,xTsfc,ro_air,Cp,
     &    undef,snow_d,vegtype)
        CALL LONGOUT(xQle,xTsfc,emiss_sfc,Stef_Boltz)
        CALL CONDUCT(icond_flag,xQc,dy_snow_z,T_old_z,gamma_z,JJ)
c        Qf = (1.0-albedo) * Qsi + Qli + xQle + xQh + xQe + xQc + Qp
        Qf = (1.0-albedo) * Qsi + Qli + xQle + xQh + xQe + xQc 
      else
        Qf = 0.0
      endif
c END J.PFLUG

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MELTTEMP(Tsfc,Tf,snow_d,vegtype)

      implicit none

      real Tsfc,snow_d,vegtype,Tf

      if (snow_d.gt.0.0 .or. vegtype.eq.20.0) then
        if (Tsfc.gt.Tf) Tsfc = Tf
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CONDUCT(icond_flag,Qc,dy_snow_z,T_old_z,gamma_z,JJ)

      implicit none

      integer icond_flag,JJ

      real Qc
      real dy_snow_z(2)
      real T_old_z(2)
      real gamma_z(2)

      if (icond_flag.eq.0) then
        Qc = 0.0
      else
        if (JJ.le.1) then
          Qc = 0.0
        else
          Qc = - (gamma_z(1) + gamma_z(2))/2.0 *
     &      (T_old_z(1) - T_old_z(2)) /
     &      (dy_snow_z(1) + dy_snow_z(2))/2.0
        endif
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CONSTS_ENBAL(Stef_Boltz,ro_air,Cp,gravity,
     &  xLs,xkappa,xLf,Tf,ro_water,Cp_water,ro_ice,emiss_sfc,
     &  Cp_snow)

      implicit none

      real Stef_Boltz,ro_air,Cp,gravity,xLs,xkappa,xLf,
     &  Tf,ro_water,Cp_water,ro_ice,emiss_sfc,Cp_snow

      emiss_sfc = 0.98
      Stef_Boltz = 5.6696e-8
      ro_air = 1.275
      Cp = 1004.
      gravity = 9.81
      xLs = 2.500e6
      xkappa = 0.4
      xLf = 3.34e5
      Tf = 273.16
      ro_water = 1000.0
      Cp_water = 4180.0
      ro_ice = 917.0
      Cp_snow = 2106.

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ENBAL(e_balance,albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,
     &  undef,snow_d,vegtype,Qp)

      implicit none

      real e_balance,albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,undef,snow_d,
     & vegtype,Qp

      if (snow_d.gt.0.0 .or. vegtype.eq.20.0) then
c        e_balance = (1.0-albedo)*Qsi+Qli+Qle+Qh+Qe+Qc-Qm+Qp 
        e_balance = (1.0-albedo)*Qsi+Qli+Qle+Qh+Qe+Qc-Qm
      else
        e_balance = undef
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE STABLEFN(stability,Tair,Tsfc,windspd,ht_windobs,
     &  gravity,xkappa,z_0)

      implicit none

      real C1,C2,B1,B2,stability,Tair,Tsfc,windspd,ht_windobs,
     &  gravity,xkappa,z_0,B8,B3,z_0_tmp

      z_0_tmp = min(0.25*ht_windobs,z_0)
      C1 = 5.3 * 9.4 * (xkappa/(log(ht_windobs/z_0_tmp)))**2 *
     &  sqrt(ht_windobs/z_0_tmp)
      C2 = gravity * ht_windobs / (Tair * windspd**2)
      B1 = 9.4 * C2
      B2 = C1 * sqrt(C2)

      if (Tsfc.gt.Tair) then
c Unstable case.
        B3 = 1.0 + B2 * sqrt(Tsfc - Tair)
        stability = 1.0 + B1 * (Tsfc - Tair) / B3
      elseif (Tsfc.lt.Tair) then
c Stable case.
        B8 = B1 / 2.0
        stability = 1.0 / ((1.0 + B8 * (Tair - Tsfc))**2)
      else
c Neutrally stable case.
        stability = 1.0
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LONGOUT(Qle,Tsfc,emiss_sfc,Stef_Boltz)

      implicit none

      real Qle,emiss_sfc,Stef_Boltz,Tsfc

      Qle = (- emiss_sfc) * Stef_Boltz * Tsfc**4

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SENSIBLE(Qh,De_h,stability,Tair,Tsfc,ro_air,Cp,
     &  undef,snow_d,vegtype)

      implicit none

      real Qh,De_h,stability,Tair,Tsfc,ro_air,Cp,undef,snow_d,vegtype

      if (snow_d.gt.0.0 .or. vegtype.eq.20.0) then
        Qh = ro_air * Cp * De_h * stability * (Tair - Tsfc)
      else
        Qh = undef
c To do this, see the laps snow-shrub parameterization runs.
c       Qh = ro_air * Cp * De_h * stability * (Tair - Tsfc)
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SFCTEMP(Tsfc,Tair,Qsi,Qli,ea,albedo,De_h,
     &  sfc_pressure,ht_windobs,windspd,ro_air,Cp,emiss_sfc,
     &  Stef_Boltz,gravity,xLs,xkappa,z_0,Tf,Qc,
     &  count_Tsfc_not_converged,prec,sprec,Cp_snow,Cp_water,xLf,Qp,
     &  dt,Twb)

      implicit none

      real Tsfc,Tair,Qsi,Qli,ea,albedo,De_h,sfc_pressure,ht_windobs,
     &  windspd,ro_air,Cp,emiss_sfc,Stef_Boltz,gravity,xLs,xkappa,
     &  z_0,Tf,Qc,AAA,CCC,DDD,EEE,FFF,C1,C2,B1,B2,z_0_tmp,
     &  count_Tsfc_not_converged,Qp,prec,sprec,Cp_snow,Cp_water,xLf,
     &  dt,Twb,rainen,snowen

c J.PFLUG
c heat advected by rain and snow precipitation is also calculated
      rainen = ((prec-sprec)*1000.0*Cp_water*max(Twb-Tsfc,0.0))/dt
      snowen = (sprec*Cp_snow*917.0*min(Twb-Tsfc,0.0))/dt
c      Qp = rainen + snowen

      AAA = ro_air * Cp * De_h
      CCC = 0.622 / sfc_pressure
      DDD = emiss_sfc * Stef_Boltz
c enable if want to include heat advected by precipitation
c      EEE = (1.0-albedo) * Qsi + Qli + Qc + Qp
      EEE = (1.0-albedo) * Qsi + Qli + Qc
      FFF = ro_air * xLs * De_h
c J.PFLUG

c Compute the constants used in the stability coefficient
c   computations.
      z_0_tmp = min(0.25*ht_windobs,z_0)
      C1 = 5.3 * 9.4 * (xkappa/(log(ht_windobs/z_0_tmp)))**2 *
     &  sqrt(ht_windobs/z_0_tmp)
      C2 = gravity * ht_windobs / (Tair * windspd**2)
      B1 = 9.4 * C2
      B2 = C1 * sqrt(C2)

      CALL SOLVE(Tsfc,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf,
     &  count_Tsfc_not_converged)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SOLVE(xnew,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf,
     &  count_Tsfc_not_converged)

      implicit none

      integer maxiter,i

      real tol,old,A,B,C,other1,other2,es0,dother1,dother2,xnew,
     &  Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf,B3,stability,
     &  dstability,fprime1,fprime2,fprime3,fprime4,B8,funct,
     &  fprime,count_Tsfc_not_converged

      tol = 1.0e-2
      maxiter = 20
      old = Tair

c Because I am interested in sublimation over snow during the
c   winter, do these calculations over ice.

c Over water.
c        A = 6.1121 * 100.0
c        B = 17.502
c        C = 240.97
c Over ice.
        A = 6.1115 * 100.0
        B = 22.452
        C = 272.55

      do i=1,maxiter

c This section accounts for an increase in turbulent fluxes
c   under unstable conditions.
        other1 = AAA * (Tair - old)
        es0 = A * exp((B * (old - Tf))/(C + (old - Tf)))
        other2 = FFF*CCC*(ea-es0)

        dother1 = - AAA
        dother2 = (- FFF)*CCC*es0*B*C/((C + (old - Tf))**2)

      if (old.gt.Tair) then
c Unstable case.
        B3 = 1.0 + B2 * sqrt(old - Tair)
        stability = 1.0 + B1 * (old - Tair) / B3
        dstability = B1/B3 - (B1*B2*(old-Tair))/
     &    (2.0*B3*B3*sqrt(old-Tair))
        fprime1 = (- 4.0)*DDD*old**3
        fprime2 = stability * dother1 + other1 * dstability
        fprime3 = stability * dother2 + other2 * dstability
        fprime4 = - 0.0

      elseif (old.lt.Tair) then
c Stable case.
        B8 = B1 / 2.0
        stability = 1.0 / ((1.0 + B8 * (Tair - old))**2)
        dstability = 2.0 * B8 / ((1.0 + B8 * (Tair - old))**3)
        fprime1 = (- 4.0)*DDD*old**3
        fprime2 = stability * dother1 + other1 * dstability
        fprime3 = stability * dother2 + other2 * dstability
        fprime4 = - 0.0

      else
c Neutrally stable case.
        stability = 1.0
        fprime1 = (- 4.0)*DDD*old**3
        fprime2 = dother1
        fprime3 = dother2
        fprime4 = - 0.0
      endif

        funct = EEE - DDD*old**4 + AAA*(Tair-old)*stability +
     &    FFF*CCC*(ea-es0)*stability +
     &    0.0
        fprime = fprime1 + fprime2 + fprime3 + fprime4

        xnew = old - funct/fprime

        if (abs(xnew - old).lt.tol) return
        old = xnew

      enddo

c If the maximum iterations are exceeded, send a message and set
c   the surface temperature to the air temperature.
c     write (*,102)
c 102 format('max iteration exceeded when solving for Tsfc')

c Count the number of times the model did not converge for this
c   time step.
      count_Tsfc_not_converged = count_Tsfc_not_converged + 1.0

      xnew = Tair

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LATENT(Qe,De_h,stability,ea,es0,ro_air,xLs,
     &  sfc_pressure,undef,snow_d,vegtype)

      implicit none

      real Qe,De_h,stability,ea,es0,ro_air,xLs,sfc_pressure,undef,
     &  snow_d,vegtype
  
      if (snow_d.gt.0.0 .or. vegtype.eq.20.0) then
        Qe = ro_air * xLs * De_h * stability *
     &    (0.622/sfc_pressure * (ea - es0))
      else
        Qe = undef
c To do this, see the laps snow-shrub parameterization runs.
c       Qe = ro_air * xLs * De_h * stability *
c    &    (0.622/sfc_pressure * (ea - es0))
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE EXCOEFS(De_h,z_0,ht_windobs,windspd,xkappa)

      implicit none

      real De_h,z_0,ht_windobs,windspd,xkappa,z_0_tmp

      z_0_tmp = min(0.25*ht_windobs,z_0)
      De_h = (xkappa**2) * windspd / ((log(ht_windobs/z_0_tmp))**2)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE VAPOR(es0,Tsfc,Tf)

      implicit none

      real es0,Tsfc,Tf,A,B,C

c Coeffs for saturation vapor pressure over water (Buck 1981).
c   Note: temperatures for Buck`s equations are in deg C, and
c   vapor pressures are in mb.  Do the adjustments so that the
c   calculations are done with temperatures in K, and vapor
c   pressures in Pa.

c Because I am interested in sublimation over snow during the
c   winter, do these calculations over ice.

c Over water.
c       A = 6.1121 * 100.0
c       B = 17.502
c       C = 240.97
c Over ice.
        A = 6.1115 * 100.0
        B = 22.452
        C = 272.55

c Compute the water vapor pressure at the surface.
      es0 = A * exp((B * (Tsfc - Tf))/(C + (Tsfc - Tf)))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GET_SFC(snow_d,albedo,z_0,Tf,Tair,snow_z0,veg_z0_tmp,
     &  vegtype,albedo_glacier,prec,albedo_diff,al_max,al_min,
     &  al_dec_cold,al_dec_melt,dt,sprec,albedo_flag)

      implicit none

      real snow_d,albedo,z_0,Tf,Tair,snow_z0,veg_z0_tmp,vegtype,
     &  albedo_veg,albedo_glacier,albedo_diff,al_max,al_min,al_dec_cold,
     &  al_dec_melt,prec,dt,tau_1,sprec

      integer albedo_flag

c The code below simulates the TIME-EVOLUTION OF SNOW ALBEDO in
c   response to precipitation and melt. 

c The general model equations are described in the paper:
c   Modeling Snow Depth for Improved Simulation of
c   Snow-Vegetion-Atmosphere Interactions.
c   Strack, J. E., G. E. Liston, and R. A. Pielke, Sr., 2004,
c   Journal of Hydrometeorology, 5, 723-734.

c J.PFLUG
c time evolution altered to mimic what had been done by 
c Sproles et al. (2013)
c If albedo is known for a POINT ONLY, can also explicitly
c define this
      if (albedo_flag.eq.0) then
        tau_1 = 86400.

c In the presence of precipitation, re-initialize the albedo.
        if (sprec.gt.0.0) then
          albedo = al_max
          if (vegtype.le.5.0) albedo = albedo - albedo_diff
        endif

c evolve the snowpack albedo
c For existing snowpack at the beginning of the timestep
        if (snow_d.gt.0.0.or.sprec.gt.0.0) then    
c define the roughness height
          z_0 = snow_z0  
c for cold snow conditions
          if (Tair.le.0.0) then
            albedo = albedo - (al_dec_cold * dt / tau_1)
            albedo = max(albedo,al_min)
c for melting snow conditions
          else 
            albedo = (albedo - al_min) * exp(-al_dec_melt * dt / tau_1)+
     &        al_min
          endif
c snow-free conditions
        else 
          if (vegtype.eq.20.0) then
            z_0 = snow_z0
            albedo = albedo_glacier
          else
            z_0 = veg_z0_tmp
            if (vegtype.le.5.0) then
              z_0 = 0.2
            elseif (vegtype.le.11.0) then
              z_0 = 0.04
            else
              z_0 = 0.02
            endif
            albedo = albedo_veg
          endif
        endif
c if point observations of albedo are known
      else
        read(98,*) albedo
        if (snow_d.eq.0.0) then
         albedo = 0.0
        endif
      endif
c END J.PFLUG
      if (albedo.lt.al_min) albedo=al_min

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE SOLVEWB_EB(xnew,Tf,Tair,rh,xLs,Cp,sfc_pressure)

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

