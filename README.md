# Introduction  
This repository is for version controlling the code that I would be writing for my undergraduate thesis. 

# Physical system of interest 
For my undergraduate thesis (with Prof. Dipankar Bhattacharya), I simulate polarization transport near neutron stars to generate pulse profiles of flux observed in desired stokes parameters. The goal behind having a model of pulse profiles is to narrow down the parameter space using observational data through bayesian inferences. Extracting a map of surface polarization using polarization observations can point us towards new developments in the multiphysics governing slab mechanisms. 

# Progress Report

Though I have officially completed the requirements of the degree, there are (primarily) two extension extensions which remains on my to-do list still <br>
1. Adding a Feautrier Module to perform radiative transfer calculations. <br>
2. Including more source functions to cover possible surface intensity distributions. <br>

May 23 - Submitted, presented and defended my thesis successfully. Thesis and bibliography available on the github repo. <br>
Apr 23 - Generated results to cover a wider region of the possible parameter space. Presented intial results at Meera Memorial Competition (at St. Stephens College) and won first prize. <br>
Mar 23 - Working on including gaussian polar cap models for surface intensity maps. Also, developing psuedo-code for inference modules and reviewing existing slab physics simulations programs. <br>
Feb 23 - Created catalog of varied NS objects, parallelized virtual 'observations', generated mock data and prepared a ten-minute summary talk (and poster). <br>
Jan 23 - Completed implementing polarization transport (as prescribed in Pavlov and Zavlin (2000) using modified Beloborodov corrections. <br>
Dec 22 - Attempted polarization transport results from Hu et al (2022). <br>
Nov 22 - Attempted analytic calculations using Sachs formalism and isentropic Schwarzschild coordinates to get the corrections. (Grad School application season = slow progress)  
Oct 22 - Gave a quarter-point update talk (chalkboard format) at Ashoka discussing elements from my thesis. Email for accessing lecture notes. <br>
Sep 22 - Understanding the theory behind how polarization undergoes parallel transport to be able to implement first and second order corrections for Schwarzschild geometries. <br>
Aug 22 - Finished the pulse profiles integration. Need to translate to the code to CUDA for GPU parallelization on university's HPC.
