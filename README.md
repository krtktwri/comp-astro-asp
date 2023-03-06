# Introduction  
This repository is for version controlling the code that I would be writing for my undergraduate thesis. 

# Physical system of interest 
For my undergraduate thesis (with Prof. Dipankar Bhattacharya), I simulate polarization transport near neutron stars to generate pulse profiles of flux observed in desired stokes parameters. The goal behind having a model of pulse profiles is to narrow down the parameter space using observational data through bayesian inferences. Extracting a map of surface polarization using polarization observations can point us towards new developments in the multiphysics governing slab mechanisms. 

# Progress Report
Now - Working on including gaussian polar cap models for surface intensity maps. Also, developing psuedo-code for inference modules and reviewing existing slab physics simulations programs. <br>
Feb 22 - Created catalog of varied NS objects, parallelized virtual 'observations', generated mock data and prepared a ten-minute summary talk (and poster). <br>
Jan 22 - Completed implementing polarization transport (as prescribed in Pavlov and Zavlin (2000) using modified Beloborodov corrections. <br>
Dec 22 - Attempted polarization transport results from Hu et al (2022). <br>
Nov 22 - Attempted analytic calculations using Sachs formalism and isentropic Schwarzschild coordinates to get the corrections. (Grad School application season = slow progress)  
Oct 22 - Gave a quarter-point update talk (chalkboard format) at Ashoka discussing elements from my thesis. Email for accessing lecture notes. <br>
Sep 22 - Understanding the theory behind how polarization undergoes parallel transport to be able to implement first and second order corrections for Schwarzschild geometries. <br>
Aug 22 - Finished the pulse profiles integration. Need to translate to the code to CUDA for GPU parallelization on university's HPC.
