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

initialize.f90: contains just one subroutine called at the beginning of the main method
that sets up the material properties set in globals.f90 and boundary/initial conditions.

GEOCHEMISTRY
=======
phreeqout/phreeqout.f90: a single box model that uses aqueous geochemistry package
PHREEQC (interfaced to fortran using IPHREEQC). presently decoupled from the fluid
dynamic model for screwing around purposes.

POST-PROCESSING
=======
*.py: scripts to generate visual output using NUMPY and MATPLOTLIB.
