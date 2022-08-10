sMass = 1.98847e30 #solar mass (kg)
c = 299792458      #speed of light (m/s)
G = 6.6743e-11     #gravitational constant (m^3/kg.s^2)

nRad = 1e4         #neutron star radius (m)
nMass = 1.6*sMass  #neutron star mass (kg)

r_G = 2*G*nMass/c**2   #Schwarzschild Radius (m)
oneover_r_G = 1/r_G    #Precomputing

lmin = 2.5981   #(3*1.73205*r_G/2 + (0.000001)) * oneover_r_G