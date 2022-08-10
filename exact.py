import numpy as np
from matplotlib import pyplot as plt, patches

#Graph Parameters
plt.rcParams['figure.figsize'] = 12,8
plt.rc('text', usetex=False)
plt.rcParams.update({'font.size': 20,
                     'legend.fontsize': 20})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['lines.linewidth'] = 2

### Coordinate converter
def X(r, theta):    
    return  r * np.cos(theta)

def Y(r, theta):
    return  r * np.sin(theta)

### Photon Trajectory Equation
def der2(r, l):          #Derivative Squared (we cannot trivially take square-root at this stage)
    d = ((r**4)/l**2)*(1 - l**2/(r**2) + r_G*l**2/r**3)
    return d

# Defining Constants
r_G = 1
lmin = 2.5981#3*np.sqrt(3)*r_G/2 + (0.000001)

### Trajectory Solver
def trajectory(l, plot=False):
        
    dphi = 0.0001   #Resolution  
    e = 0.000001    #Tolerance
    
    sign = -1
    zero = False 
    
    # Input Array
    phiArray = np.arange(0, 5*np.pi, dphi)
    
    rArray = []   #List for storing solution (will later turn into numpy array for easier manipulations)
    rArray.append(1000)  #Initial condition (Ray starts virtually at infinity)
    
    
    #### Actual trajectory computation begins here ####
    i = 1
    while i<len(phiArray):
        if zero == False:
            if der2(rArray[i-1], l) < e:
                zero = True
                sign = +1
            else:
                zero = False
        
        if rArray[i-1]<20000 and rArray[i-1]>r_G:
            
            #Runge-Kutta 4th Order
            k1 = dphi * np.sqrt(der2(rArray[i-1], l))
            k2 = dphi * np.sqrt(der2(rArray[i-1]+ k1/2, l))
            k3 = dphi * np.sqrt(der2(rArray[i-1]+ k2/2, l))
            k4 = dphi * np.sqrt(der2(rArray[i-1]+ k3, l))
            
            rArray.append(rArray[i-1] + sign * (k1 + 2*k2 + 2*k3 + k4)/6)
        else:
            break
        
        i+=1
    #### Trajectory computation ends here ####
    
    #Polar coordinates to cartersian coordinates
    x = X(rArray[0:i], phiArray[0:i])
    y = Y(rArray[0:i], phiArray[0:i])
    
    if plot==True:
        fig = plt.figure()
        ax = fig.add_subplot()
        plt.gca().set_aspect('equal')

        
        ax.plot(x, y)
        
        circle = patches.Circle((0, 0), radius=r_G, color='black', zorder=100)
        ax.add_patch(circle)
        ax.patch.set_visible(False)
        ax.set_xlabel(r"X/$r_G$")
        ax.set_ylabel(r"Y/$r_G$")

        plt.xlim([-10, 10])
        plt.ylim([-10, 10])
    
    return x, y




