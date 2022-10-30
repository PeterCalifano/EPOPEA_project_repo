SPICE-Documentation:
	-SPICEforMATLABHAndbook: guide for user to install SPICE tool on Matlab 	and get familiar with it

	-SPICE_overview and SPICE_overview_Addendum: slides provided by GNC 		course. There are some of the main commands with examples and exercises

	- see NASA site: https://naif.jpl.nasa.gov/naif/tutorials.html for 	 	tutorials

KERNELS:	
	-de440s.bsp: ephemerides of planets wrt SSB. coverage 1849-2150
	-sat441.bsp: ephemerides of moons of Saturn wrt Saturn Bar. coverage 				 1749 - 2250
	-naif0012.tls: leap second kernel
	-pck00010_mod.tpc: physical kernel containing data of several objects
		             !! Modification: Enceladus radii are all equal!!
	-LATs.def: definition kernel for creation of ephemerides kernel LATs.bsp 		     and frame kernel LATs.tf