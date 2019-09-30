# OPF-based-Under-Frequency-Load-Shedding

## Data files:
	* busdata.txt: file containing the bus data of the IEEE 5-bus power system;
	* linedata.txt: file containing the line data of the IEEE 5-bus power system

## Scripts and functions needed to solve the UFLS optimzation problem:
	* startup_RTS: called by the initialization function in Simulink. Set the base values, compute the lines lengths and set the load power profile and the amount of load shedding in each load;
	* startup_UFLS.m: call 'createGridFromFiles.m' and generate the 'UFLS' structure; 
	* createGridFromFiles.m: extract data from 'busdata.txt' and 'linedata.txt' and call 'YMatrix.m';
	* YMatrix.m: compute the admittance matrix from 'linedata.txt';
	* UFLS_optimization_problem: call 'startup_UFLS.m' to get the 'UFLS' structure and call 'sensCoeffs.m' to compute the sensitivity coefficients. Then, solve the optimization problem and save the results;
	* sensCoeffs.m: compute the sensitivity coefficients

## Simulink models:
	* IEEE5BusSystemModel_opt_UFLS.slx: Simulink model of the IEEE 5-bus power system with optimized UFLS;
	* IEEE5BusSystemModel_std_UFLS.slx: Simulink model of the IEEE 5-bus power system with standard UFLS scheme

## Miscellaneous:
	* setup_fitting.m: setup of time and frequency vectors in order to use Matlab fitting tools to fit the SFR model parameters;
	* SFRmodel.sfit: fitting session for the SFR model;  
	* sensitivity_coefficients_derivation.docx: derivation of the sensitivity coefficients considering voltage-dependent loads
