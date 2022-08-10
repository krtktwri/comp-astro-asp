import numpy as np
from modules.myconstants import r_G

### Coordinate converter
def X(r, theta):    
    return  r * np.cos(theta)

def Y(r, theta):
    return  r * np.sin(theta)

def r(x, y):
    return np.sqrt(x**2 + y**2)

def theta(x, y):
    return np.arctan(y/x)

'''
Important: For the following functions, both input and output are normalized wrt r_G
'''
### Conversion between angle ray angle and impact parameter (beloborodov2002 Eq. A2)
def l2alpha(l, R):
    l = l*r_G
    R = R*r_G
    return np.arcsin((l/R)*np.sqrt(1 - r_G/R))
    
### Conversion between angle ray angle and impact parameter (beloborodov2002 Eq. A2)
def alpha2l(alpha, R):
    R = R*r_G
    return ((np.sin(alpha)*R)/np.sqrt(1 - r_G/R))/r_G