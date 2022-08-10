import numpy as np
from matplotlib import pyplot as plt, patches

#my packages
from modules.myconstants import *
from modules.converter import *

#Graph Parameters
plt.rcParams['figure.figsize'] = 12,8
plt.rc('text', usetex=False)
plt.rcParams.update({'font.size': 15,
                     'legend.fontsize': 15})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 15
plt.rcParams['lines.linewidth'] = 2

### Photon Trajectory Equation (derivation - bhattacharya2003)

def der2(r, l):          #Derivative Squared (we cannot trivially take square-root at this stage)
    d = ((r**4)/l**2)*(1 - l**2/(r**2) + r_G*(l**2)/r**3)
    return d

### Trajectory Solver
def trajectory(l, rad = r_G, plot=False):
    
    # Photon orbits circle around the black hole many times. 2pi rotation would not complete the trajectory 
    if l == lmin:
        factor = 8
    if l != lmin:
        factor = 2
    
    dphi = 1e-4    #Resolution  
    e = l      #Tolerance
    
    # Input Array
    phiArray = np.arange(0, factor*np.pi, dphi)
    
    l = l*r_G
    
    stop = rad

    sign = -1
    zero = False 
    
    
    rArray = []   #List for storing solution (will later turn into numpy array for easier manipulations)
    rArray.append(1000*r_G)  #Initial condition (Ray starts virtually at infinity) 
    
    
    #### Actual trajectory computation begins here ####
    i = 1
    while i<len(phiArray):

        if zero == False:
            if der2(rArray[i-1], l) < e:
                zero = True
                sign = +1
            else:
                zero = False         
         
        if rArray[i-1]<1100*r_G and rArray[i-1]>=stop:
            
            #Runge-Kutta 4th Order
            k1 = dphi * np.sqrt(der2(rArray[i-1], l))
            k2 = dphi * np.sqrt(der2(rArray[i-1]+ k1/2, l))
            k3 = dphi * np.sqrt(der2(rArray[i-1]+ k2/2, l))
            k4 = dphi * np.sqrt(der2(rArray[i-1]+ k3, l))
            
            rArray.append((rArray[i-1] + sign * (k1 + 2*k2 + 2*k3 + k4)/6))
        else:
            break
            
        i+=1       
    #### Trajectory computation ends here ####
    
    #Polar coordinates to cartersian coordinates
    x = X(rArray[0:i], phiArray[0:i])
    y = Y(rArray[0:i], phiArray[0:i])
    
    ### Optional Quick Plotting 
    if plot==True:
        fig = plt.figure()
        ax = fig.add_subplot()
        plt.gca().set_aspect('equal')

        
        ax.plot(x*oneover_r_G, y*oneover_r_G, label=f'Impact Param {np.round(l/r_G, 4)}')
        
        surface = patches.Circle((0, 0), radius=stop*oneover_r_G, label="Surface", color='blue', alpha=0.3, zorder=100)             
        horizon = patches.Circle((0, 0), radius=r_G*oneover_r_G, label="Horizon", color='black', alpha=0.3, zorder=200)

        ax.add_patch(surface)
        ax.add_patch(horizon)
        ax.patch.set_visible(False)
        ax.set_xlabel(r"X/$r_G$")
        ax.set_ylabel(r"Y/$r_G$")
        
        plt.legend()
        
        lim = 10

        plt.xlim([-lim, lim])
        plt.ylim([-lim, lim])
        
    
    return x*oneover_r_G, y*oneover_r_G




