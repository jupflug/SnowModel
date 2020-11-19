# SnowModel development: Gravity-drainage snowpack liquid water percolation

![combined](https://user-images.githubusercontent.com/20308427/92159459-4a1c9380-ede2-11ea-8a07-99015fa9e4a9.jpg)

##### This repository contains developments to the SnowModel distributed model framework [(Liston and Elder, 2004)](https://journals.ametsoc.org/jhm/article/7/6/1259/5465) detailed in [Pflug et al. (2019)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018WR024632@10.1002/(ISSN)1944-7973.SNOWEX1). Model developments were focused on liquid water movement in the modeled snowpack and are appropriate for only multilayer simulations. Liquid water movement relies more-heavily on snow temperature than the default model and therefore requires simulations with model forcing at at least 3-hourly timesteps, or more frequent. In addition to developments for liquid water percolation:
* An option for the user to explicitly define frozen precipitation (instead of using the default precipitation-phase partition) is included
* The 2-layer model used in Pflug et al. (2019) is included 

### Users are cautioned that the developments here are not yet included in the publicly-available SnowModel source code. Differences in model versions are possible. Users are encouraged to inquire about the most-recent model developments and whether the model here is best for their purposes.

## Getting started
* Download or clone this repository
* Move all files except snowmodel.par in a directory named /code
* Organize forcing, outputs, and supplementary data as directed by snowmodel.par
* Use the compiler (pgf77, gfortran, g77, f77, etc.) of choice to compile and run SnowModel 

###### Sources
* Liston, G.E. and Elder, K.,2006. A Distributed Snow-Evolution Modeling System (SnowModel). J. Hydrometeor. 7, 1259-1276.
* Pflug, J.M., Liston, G.E., Nijssen, B., Lundquist J.D., 2019. Testing Model Representations of Snowpack Liquid Water Percolation Across Multiple Climates. Water Resources Research. 55(6), 4820-4838.

