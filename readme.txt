This directory contains the MATLAB code for solving the coastal-trapped wave problem from `Spectral methods for coastal-trapped waves and instabilities in a background flow' by by M. N. Crowe and E. R. Johnson.

Scripts for the three examples considered in this paper are provided. These files have been tested on Matlab 2022b-2023b and are expected to work on  some older versions. All scripts are written by M. N. Crowe based on theory devised by M. N. Crowe and E. R. Johnson. For comments, questions or suggestions, please email Matthew.Crowe2@newcastle.ac.uk.

List of files:

functions/				- Directory containing all functions required to solve the CTW problem
	cmap2.m				- Creates custom red/white/blue colorbar, dependency of Plot_Mode.m
	Create_Grid.m			- Creates the grid structure input for Find_Modes.m
	create_operator.m		- Builds a multi-dimensional spectral operator, dependency of Find_Modes.m
	Create_Params.m			- Creates the parameters input for Find_Modes.m
	Find_Modes.m			- Finds the frequency/growth rate and eigenfunctions for a given CTW/instability problem
	Gen_EVP.m			- Solves a generalised eigenvalue problem using MATLAB `eigs', dependency of Find_Modes.m
	grid_composite.m		- Creates a composite grid, dependency of Create_Grid.m
	grid_spectral.m			- Creates a spectral grid and differentiation operator, dependency of Create_Grid.m
	Plot_Mode.m			- Creates a color plot of a 2D function, used for plotting examples only
Barotropic_Instability/			- Directory containing all functions required to solve the barotropic problem (Example 2)
	Create_Grid_Barotropic.m	- Creates the grid input for Find_Modes_Barotropic.m and Find_Modes_Equiv_Barotropic.m
	Create_Params_Barotropic.m	- Creates the parameters input for Find_Modes_Barotropic.m and Find_Modes_Equiv_Barotropic.m
	Find_Modes_Barotropic.m		- Finds the frequency/growth rate and eigenfunctions for the barotropic problem
	Find_Modes_Equiv_Barotropic.m	- Finds the frequency/growth rate and eigenfunctions for the equivalent barotropic problem
Example_1_Kelvin_Waves_on_Shelf.m	- Idealised Kelvin and shelf wave example from Sec. 4.1
Example_2_Quasi_Barotropic_Jet.m	- Idealised jet instability example from Sec. 4.2
Example_3_Realistic_CTW_Modes.m		- Realistic CTWs example from Sec. 4.3
Gelderloos_et_al_dat.mat		- Data file containing bottom topography profile from Gelderloos et. al. (2021)




