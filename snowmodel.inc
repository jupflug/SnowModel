c snowmodel.inc

      integer nx_max,ny_max,nstns_max,nvegtypes,max_time_steps,
     &  max_obs_dates

      parameter (nx_max=1501,ny_max=1501)
      parameter (nstns_max=10000)
      parameter (nvegtypes=30)
      parameter (max_time_steps=130000)

c max_obs_dates is used in the data assimilation routines.  It's
c   value must be greater than the number of observation dates
c   in the entire simulation + (plus) the number of years in the
c   simulation.  For example, for a 6-year simulation with two
c   observation dates in each year, you would set max_obs_dates
c   to be at least = 18.
      parameter (max_obs_dates=12)

c If you are running the multi-layer snow model (even with a single
c   layer) nz_max must be at least one greater than max_layers in
c   snowmodel.par.  This is because the model will build a new layer
c   with the new snowfall and then it is merged with the layer below
c   if you only want a single snow layer.  If you are running
c   SnowModel's original single layer model, nz_max can be 1 (but if
c   nz_max=2 it will avoid a warning message if you are compiling
c   the code with gfortran).
      integer nz_max

      parameter (nz_max=25)
c     parameter (nz_max=2)

