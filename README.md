MASSACR
=======
Modeling Altered Seafloor: Simulation and Climatic Response

BACKGROUND
=======
MASSACR improves upon its predecessor, SWAGCUNT, with more adaptable PDE solvers, cleaner
output, and improved functionality for all of your porous-media-flow-coupled-to-aqueous-geochemistry 
needs. MASSACR incorporates new insights (anderson et al. 2012, 2013) into off-axis hydrothermal systems
that have not been previously considered in a geochemical context. not that i know of anyway.

FLUID DYNAMICS
=======
massacre.f90: main method, PDE solving subprograms, message passing central hub (distribution to slave
processors and collection by root processor), output writing. reactive transport happens here (done
by root processor) instead of in its own module. 

globals.f90: includes stuff that the model needs to be available globally, including meshgrid 
points and timesteps in various fine+coarse combinations, physical properties of fluid/rock, and solvers 
for sets of linear equations (e.g. banded matrix solver, tridiagonal solver)

initialize.f90: contains just one subroutine, init(), called at the beginning of the 
main method that sets up the material properties set in globals.f90 and boundary/
initial conditions. fluid/rock properties that require extensive calculation are born here
instead of the globals module.

GEOCHEMISTRY
=======

alteration.f90: this is a module that contains all of the necessary functionality for
basalt alteration that can be coupled to the fluid dynamic model by running alter()
for one or all or some grid cell(s). it is designed to accept relevant information 
such as cell temperature, timestep size, amounts of primary consituents, and amounts of
secondary alteration products which are then used as inputs to PHREEQC. the output is
PHREEQC SELECTED-OUTPUT and is not parsed in alter().

POST-PROCESSING
=======

plotPressure.py: plot streamfunctions, isotherms, and several benchmark values that 
correspond to output in hydrothermal circulation simulations by snelgrove and forster,
1996.

plotAnimation.py: for making animations for talks and navahnavahnavah.com

ONE-BOX MODELS
=======

i made a bunch of models for screwing around purposes. each does something slightly
different and the directory contains individual post-processing and visualization 
scripts.

phreeqout/: original basalt alteration experiments (based on carbon sequestration 
simulations by pham et al., 2012) using aqueous geochemistry package PHREEQC 
interfaced to FORTRAN using IPHREEQC.

batch/: single box model of basalt alterating and carbonate precipitation including
mixing of seawater-derived hydrothermal fluid with fresh deep ocean seawater.
batch.f90 and batchControl.f90 handle the mixing differently. batchControl.f90 is
probably better.

cell/: pretty much the same as batchControl but with better output/post-processing
and ultimately for the development of MASSACR's geochemical model.
