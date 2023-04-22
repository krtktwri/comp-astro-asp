import numpy as np

# Defining a class for details of the polar cap
class capDetails:
    def __init__(self, r, NM, base=0):
        self.r = r         # angular span of the polar cap (in degrees)
        self.I_o = NM[0]   # intensity of the polar cap in Ordinary Mode
        self.I_e = NM[1]   # intensity of the polar cap in Extraordinary Mode
        self.base = base

# Defining an intensity function 
def gaussianCap(phi, theta, r, base, eta):

    # note the polar cap center always lies on the X-Z plane
    # converting polar cap angular size to radians
    r = np.deg2rad(r) 
        
    # gaussian intensity distribution at the poles
    if (phi >= (eta-r) and phi <= (eta+r)):
        if np.abs(theta) <= r:
            return np.exp(-(phi-eta)**2/(2*r**2)) * np.exp(-(theta)**2/(2*r**2))
        else:
            return base
        
    # diametrically opposing hotspot
    if (phi >= (eta+np.pi-r) and phi <= (eta+np.pi+r)):
        if np.abs(theta) <= r:
            return np.exp(-(phi-(eta+np.pi))**2/(2*r**2)) * np.exp(-(theta)**2/(2*r**2))
        else:
            return base
    
    # low intensity at the equator
    else:   
        return base