# Introduction  
This repository is for version controlling the code that I am writing for my undergraduate thesis. 

# Physical system of interest 
For my undergraduate thesis, with my advisor (Prof. Dipankar Bhattacharya), I simulate polarization transport around neutron stars to generate pulse profiles of stokes parameters. These pulse profiles are fed into bayesian inference programs to extract a surface polarization map from the observational data. Understanding the polarization on the surface of a neutron star can give us insights about the slab physics governing these ultra-dense objects and point towards interesting developments in QED.

# Updates
Now - On computational side, exploring GPU based inference routines. On theoretical side, studying about vacuum birefringence and slab physics. 
Jan 22 - Finished implementing polarization transport as described in Pavlov and Zavlin (2000). Can produce stokes parameter intensity profiles now.
Dec 22 - Understanding polarization transport results from Hu et al (2022).
Nov 22 - Attempted analytic calculations using Sachs formalism and isentropic Schwarzschild coordinates to get the corrections. (Grad School application season = slow progress)
Oct 22 - Gave a quarter-point update talk (chalkboard format) at Ashoka discussing elements from my thesis. Email for accessing lecture notes.
Sep 22 - Understanding the theory behind how polarization undergoes parallel transport to be able to implement first and second order corrections for Schwarzschild geometries.
Aug 22 - Finished the pulse profiles integration. Need to translate to the code to CUDA for GPU parallelization on university's HPC.