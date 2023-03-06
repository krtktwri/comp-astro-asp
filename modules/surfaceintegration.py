import numpy as np
from modules.approximations import *
from modules.myconstants import *


'''
Geometry adapted from Appendix 2 of mb2011 (see '/papers' folders)
'''

def psi(theta, phi, c, phase): 
    
    i = c[1]
    w = phase
        
    n_psi_hat = (np.sin(i)*np.sin(w), np.sin(i)*np.cos(w), np.cos(i))                               #Unit vector towards lOS 
    r_hat = (np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)) #Unit vector to a position on surface 
    
    return np.arccos(np.dot(n_psi_hat, r_hat)) #Angle between n_psi_hat and r_hat

# Normalized Intensity Distribution
def I(theta, phi, c):
    
    eta_p = c[0]
    theta_p = np.pi/2
    r = 0.1 #radius in radians
    
    # Primary Pole
    if (np.abs(phi-eta_p)<= r and np.abs(theta-theta_p)<= r):
        return 5
    
    #Secondary Pole
    elif (np.abs(phi-(eta_p + np.pi))<= r and np.abs(theta-(theta_p))<= r):
        return 3
    elif (np.abs(phi-(eta_p))<= r and np.abs(theta-(theta_p+np.pi))<= r):
        return 3
    
    else:
        return 0.1
    
# Flux per unit area (Eq. 18 - mb2011)
def flux_per_unitArea(theta, phi, c, phase, dA):
    
    a = psi2alphaCorrected(psi(theta, phi, c, phase), R)   
    
    if a>=np.pi/2:
        return 0
    else:
        return I(theta, phi, c) * (dA) * np.cos(a)
    


# Acquiring complete pulse profile    
def pulseProfile(phase, c, res='low'): 
    
    if res == 'low':
        dtheta, dphi = 0.1, 0.1        # Low Resolution for testing
    else:
        dtheta, dphi = 0.05, 0.05      # Higher Resolution for final
    
    # Initializing Arrays
    thetaRange = np.arange(-np.pi, np.pi, dtheta)
    phiRange = np.arange(0, 2*np.pi, dphi)

    F = []

    #Summing over all area patches on the surface
    for i in range(0, len(thetaRange)):
        for j in range(0, len(phiRange)):
            F.append(flux_per_unitArea(thetaRange[i], phiRange[j], c, phase, (R**2*dtheta*dphi)))
        
    return sum(F)            