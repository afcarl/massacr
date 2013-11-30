MASSACR
=======
Modeling Altered Seafloor: Simulation and Climatic Response

BACKGROUND
=======
MASSACR improves upon its predecessor, SWAGCUNT, with cleaner output, faster solvers, 
and improved functionality for all of your porous-media-flow-coupled-to-aqueous-
geochemistry needs.

FLUID DYNAMICS
=======
massacre.f90: includes the main method and the subprograms that solve the PDEs that
describe porous media flow

globals.f90: includes specifications that aren't subject to change later 
(e.g. grid resolution, timestep size, physical properties of fluid) and solvers 
for sets of linear equations (e.g. banded matrix solver, tridiagonal solver)

initialize.f90: contains just one subroutine, init(), called at the beginning of the 
main method that sets up the material properties set in globals.f90 and boundary/
initial conditions.

GEOCHEMISTRY
=======
phreeqout/: this is where i originally developed the basalt alteration experiments
(based on carbon sequestration simulations by pham et al., 2012) by interfacing
aqueous geochemistry package PHREEQC to FORTRAN using IPHREEQC. everything in
this directory is decoupled from the fluid dynamic model for screwing around
purposes.

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

phreeqout/plotGlass.py: parse and visualize basalt alteration, pH evolution, secondary 
mineral precipitation and carbonate precipitation in carbon sequestration experiments 
based off those by pham et al., 2012.
