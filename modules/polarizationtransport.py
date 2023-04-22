import numpy as np 
from modules.approximations import *

from modules.surfaceintegration import *
from modules.myconstants import *
from modules.surfacemap import *


'''
Function to calculate degree of linear polarization
Input: Angle between the magnetic field and the scattering plane (theta_prime), magnetic field strength (B), and energy of photon (E)
Output: Degree of linear polarization 
'''
def linear_polarization(theta_prime, B, E):
    
    E_Be = 11.6 * B # Electron cyclotron energy in kev (with B expressed in 10^12 gauss)
    E_Bi = 0.00632 * B # Ion cyclotron energy in kev (with B expressed in 10^12 gauss)
     
    q = (E**2 *(E_Be**2 + E_Bi**2 - E_Be*E_Bi) - E_Be**2*E_Bi**2)/(E**3*(E_Be - E_Bi)) 
    
    return np.abs(q)*np.sin(theta_prime)**2/ np.sqrt(4*np.cos(theta_prime)**2 + q**2*np.sin(theta_prime)**4)

'''
Function to calculate the total flux of a star in a given phase
Input: Phase, star object, energy, resolution
Output: Total flux of the star in that phase
'''
def F_I(phase, star, E, res='low'):
    
    I_o, I_e = 1, 0
    
    # unpacking star parameters
    R, star_eta, star_i = star.radius, star.eta, star.i
    
    if res == 'low':
        dtheta, dphi = 0.1, 0.1        # Low Resolution for testing
    else:
        dtheta, dphi = 0.05, 0.05      # Higher Resolution for final

    phiLim = alpha2psiCorrected(np.pi/2, R)
    # Initializing Arrays
    thetaRange = np.arange(0, np.pi, dtheta)
    phiRange = np.arange(-phiLim, phiLim, dphi)

    dA = R**2*dtheta*dphi       # unit area element
    u = 1/R                    # ratio R_g/R (R already in terms of schwarzschild radius)   
    g_R = np.sqrt(1 - (u))     # gravitational redshift parameter
    
    F = []                     # empty list to store flux for each patch
    
    # from each alpha and theta
    for i in range(0, len(thetaRange)):
        for j in range(0, len(phiRange)):
            
            ## calculate psi (using belo function)
            a = psi2alphaCorrected((np.pi/2 - phiRange[j]), R)         
                        
            ## compute flux in Q stokes parameter from this area element
            F.append(g_R*dA*np.cos(a)*(I_o+I_e))
            
    return np.sum(F)

'''
Function to integrate stokes parameter over the visible part of the star (with even constant intensity distribution)
Input: phase, star object, energy, resolution (default = 'low')
Ouput: FQ, FU (stokes parameters)
'''
def stokes(phase, star, E, res='low'): 
    
    I_o, I_e = 1, 0
    
    if res == 'low':
        dtheta, dphi = 0.1, 0.1        # Low Resolution for testing
    else:
        dtheta, dphi = 0.05, 0.05      # Higher Resolution for final
        
    # unpacking star parameters
    R, star_eta, star_i, B_p = star.radius, star.eta, star.i, star.B_p     
    
    # extracting correction coefficients
    corr_coeff = extract_coefficients(R) 
    phiLim = np.deg2rad(160)   # even in the most compact neutron stars (R = 2), lensing doesn't lead to seeing patches beyond 160deg
    
    # Initializing Arrays
    thetaRange = np.arange(0, np.pi, dtheta)
    phiRange = np.arange(-phiLim, phiLim, dphi)

    dA = R**2*dtheta*dphi       # unit area element
    Theta = np.arccos(np.cos(star_eta)*np.cos(star_i) + np.sin(star_eta)*np.sin(star_i)*np.cos(phase))    # angle between LOS and magnetic axis [WARNING: Theta is not the same as theta AND I am not sure whether this is the correct relation]
    
    u = 1/R                    # ratio R_g/R (R already in terms of schwarzschild radius)   
    g_R = np.sqrt(1 - (u))     # gravitational redshift parameter
        
    # GR correction to magnetic dipole moment (eq. 14 - pavlov2000)
    f = 2*(u**2 - 2*u - 2*(1-u)*np.log(1-u))/(np.sqrt(1-u)*(u**2 + 2*u + 2*np.log(1-u)))
    
    FQ, FU = [], []                     # empty list to store flux for each patch
    
    # from each alpha and theta
    for i in range(0, len(thetaRange)):
        for j in range(0, len(phiRange)):
            
            ## calculate psi (using belo function)
            a = fastApprox((np.pi/2 - phiRange[j]), R, corr_coeff)[0]         
            
            if a>=np.pi/2:
                pass
            
            ### r_hat vector components in x, y, z coordinates (eq 15 - pavlov2000)
            r_hat = np.array([np.sin(phiRange[j])*np.cos(thetaRange[i]), np.sin(phiRange[j])*np.sin(thetaRange[i]), np.cos(phiRange[j])])
                        
            ### Defining component transformation matrix 
            T = np.array([[np.cos(thetaRange[i])*np.cos(phiRange[j] - a), np.sin(thetaRange[i])*np.cos(phiRange[j] - a), -np.sin(phiRange[j] - a)],
                           [-np.sin(thetaRange[i]), np.cos(thetaRange[i]), 0],
                           [np.cos(thetaRange[i])*np.sin(phiRange[j] - a), np.sin(thetaRange[i])*np.sin(phiRange[j] - a), np.cos(phiRange[j] - a)]])
            
            ### r_hat vector components  in x', y', z' coordinates (eq 8-10 - pavlov2000)
            # r_hat_prime = np.dot(T, r_hat) 
            
            ## Defining m_hat vector components in x, y, z coordinates (eq 16 - pavlov2000)
            m_hat = np.array([np.sin(Theta), 0, np.cos(Theta)])
            
            ## get magnetic field at all points (eq 13 - pavlov2000) 
            B_vec = (B_p/2)*((2+f)*np.dot(r_hat, m_hat)*r_hat - f*m_hat)  # in the basis vectors of (x,y,z)
            B_vec_prime = np.dot(T, B_vec)                                # in the prime coordinates (x', y', z')
            
            ## find phi' using By' and Bx'.  
            phi_prime = np.arctan2(B_vec_prime[1], B_vec_prime[0])               # phi' is the angle between the x' axis and the projection of B_vec_prime onto the x'y' plane

            ## angle between magnetic field and the unit wave vector at the surface
            Theta_prime = np.arccos(B_vec_prime[2]/np.linalg.norm(B_vec_prime))  
            
            ## angle between B_vec and normal to the surface
            Theta_B = np.arccos(np.dot(B_vec/np.linalg.norm(B_vec), r_hat))
            B = B_p * f / np.sqrt(4 - (4 - f**2)*np.cos(Theta_B)**2)  # magnetic field strength at the surface (eq. 17 - pavlov2000)
            
            # degree of linear polarization
            p_L = linear_polarization(Theta_prime, B, g_R*E)
            
            ## compute flux in Q stokes parameter from this area element
            FQ.append(g_R*dA*np.cos(a)*(I_o-I_e)*p_L*np.cos(2*(thetaRange[i] - phi_prime)))
            FU.append(g_R*dA*np.cos(a)*(I_e-I_o)*p_L*np.sin(2*(thetaRange[i] - phi_prime)))

    return np.nansum(FQ), np.nansum(FU)


###############################################################
#### INTRODUCING GAUSSIAN POLAR CAP SURFACE MAP ###############
###############################################################

def F_I_PolarCap(phase, star, E, cap, res='low'):
        
    # unpacking star parameters
    R, star_eta, star_i = star.radius, star.eta, star.i
    r, I_o, I_e, base = cap.r, cap.I_o, cap.I_e, cap.base 

    if res == 'low':
        dtheta, dphi = 0.1, 0.1        # Low Resolution for testing
    else:
        dtheta, dphi = 0.05, 0.05      # Higher Resolution for final

    phiLim = alpha2psiCorrected(np.pi/2, R)
    
    # Initializing Arrays
    thetaRange = np.arange(0, np.pi, dtheta)
    phiRange = np.arange(-phiLim, phiLim, dphi)

    dA = R**2*dtheta*dphi       # unit area element
    u = 1/R                    # ratio R_g/R (R already in terms of schwarzschild radius)   
    g_R = np.sqrt(1 - (u))     # gravitational redshift parameter
    
    F = []                     # empty list to store flux for each patch
    
    # from each alpha and theta
    for i in range(0, len(thetaRange)):
        for j in range(0, len(phiRange)):
            
            I_o, I_e = gaussianCap(phiRange[j], thetaRange[i], r, base, star_eta) *np.array([I_o, I_e])

            ## calculate psi (using belo function)
            a = psi2alphaCorrected((np.pi/2 - phiRange[j]), R)         
                        
            ## compute flux in Q stokes parameter from this area element
            F.append(g_R*dA*np.cos(a)*(I_o+I_e))
            
    return np.sum(F)


'''
Stokes parameter calculation but with gaussian polar cap surface map
'''
def stokes_PolarCap(phase, star, E, cap, res='low'): 
        
    if res == 'low':
        dtheta, dphi = 0.1, 0.1        # Low Resolution for testing
    else:
        dtheta, dphi = 0.05, 0.05      # Higher Resolution for final
        
    # unpacking star parameters
    R, star_eta, star_i, B_p = star.radius, star.eta, star.i, star.B_p    
    r, I_o, I_e, base = cap.r, cap.I_o, cap.I_e, cap.base 
    
    phiLim = alpha2psiCorrected(np.pi/2, R)
    # Initializing Arrays
    thetaRange = np.arange(0, np.pi, dtheta)
    phiRange = np.arange(-phiLim, phiLim, dphi)

    dA = R**2*dtheta*dphi       # unit area element
    Theta = np.arccos(np.cos(star_eta)*np.cos(star_i) + np.sin(star_eta)*np.sin(star_i)*np.cos(phase))    # angle between LOS and magnetic axis [WARNING: Theta is not the same as theta AND I am not sure whether this is the correct relation]
    
    u = 1/R                    # ratio R_g/R (R already in terms of schwarzschild radius)   
    g_R = np.sqrt(1 - (u))     # gravitational redshift parameter
        
    # GR correction to magnetic dipole moment (eq. 14 - pavlov2000)
    f = 2*(u**2 - 2*u - 2*(1-u)*np.log(1-u))/(np.sqrt(1-u)*(u**2 + 2*u + 2*np.log(1-u)))
    
    FQ, FU = [], []                     # empty list to store flux for each patch
    
    # from each alpha and theta
    for i in range(0, len(thetaRange)):
        for j in range(0, len(phiRange)):
            
            I_o, I_e = gaussianCap(phiRange[j], thetaRange[i], r, base, star_eta) * np.array([I_o, I_e])
            
            ## calculate psi (using belo function)
            a = psi2alphaCorrected((np.pi/2 - phiRange[j]), R)         
        
            ### r_hat vector components in x, y, z coordinates (eq 15 - pavlov2000)
            r_hat = np.array([np.sin(phiRange[j])*np.cos(thetaRange[i]), np.sin(phiRange[j])*np.sin(thetaRange[i]), np.cos(phiRange[j])])
                        
            ### Defining component transformation matrix 
            T = np.array([[np.cos(thetaRange[i])*np.cos(phiRange[j] - a), np.sin(thetaRange[i])*np.cos(phiRange[j] - a), -np.sin(phiRange[j] - a)],
                           [-np.sin(thetaRange[i]), np.cos(thetaRange[i]), 0],
                           [np.cos(thetaRange[i])*np.sin(phiRange[j] - a), np.sin(thetaRange[i])*np.sin(phiRange[j] - a), np.cos(phiRange[j] - a)]])
            
            ### r_hat vector components  in x', y', z' coordinates (eq 8-10 - pavlov2000)
            # r_hat_prime = np.dot(T, r_hat) 
            
            ## Defining m_hat vector components in x, y, z coordinates (eq 16 - pavlov2000)
            m_hat = np.array([np.sin(Theta), 0, np.cos(Theta)])
            
            ## get magnetic field at all points (eq 13 - pavlov2000) 
            B_vec = (B_p/2)*((2+f)*np.dot(r_hat, m_hat)*r_hat - f*m_hat)  # in the basis vectors of (x,y,z)
            B_vec_prime = np.dot(T, B_vec)                                # in the prime coordinates (x', y', z')
            
            ## find phi' using By' and Bx'.  
            phi_prime = np.arctan2(B_vec_prime[1], B_vec_prime[0])               # phi' is the angle between the x' axis and the projection of B_vec_prime onto the x'y' plane

            ## angle between magnetic field and the unit wave vector at the surface
            Theta_prime = np.arccos(B_vec_prime[2]/np.linalg.norm(B_vec_prime))  
            
            ## angle between B_vec and normal to the surface
            Theta_B = np.arccos(np.dot(B_vec/np.linalg.norm(B_vec), r_hat))
            B = B_p * f / np.sqrt(4 - (4 - f**2)*np.cos(Theta_B)**2)  # magnetic field strength at the surface (eq. 17 - pavlov2000)
            
            # degree of linear polarization
            p_L = linear_polarization(Theta_prime, B, g_R*E)
            
            ## compute flux in Q stokes parameter from this area element
            FQ.append(g_R*dA*np.cos(a)*(I_o-I_e)*p_L*np.cos(2*(thetaRange[i] - phi_prime)))
            FU.append(g_R*dA*np.cos(a)*(I_e-I_o)*p_L*np.sin(2*(thetaRange[i] - phi_prime)))

    return np.nansum(FQ), np.nansum(FU)

# create a class for observation
'''
Am observation class represents a single observation of a given neutron star, 
for a specific phase range (in radians), at a given energy.
'''
class observation:
    def __init__(self, star, phaseRange, E, I_map=[]):
        
        if I_map == []:
            stokesParamComputed = [stokes(m, star, E, res='low') for m in phaseRange]
            stokesOutput = np.transpose(stokesParamComputed)
            totalI = [F_I(m, star, E, res='low') for m in phaseRange]
        else:
            stokesParamComputed = [stokes_PolarCap(m, star, E, I_map, res='low') for m in phaseRange]
            stokesOutput = np.transpose(stokesParamComputed)     
            totalI = [F_I_PolarCap(m, star, E, I_map, res='low') for m in phaseRange]
       
        # Computed Stokes Parameters
        self.FQ = stokesOutput[0]
        self.FU = stokesOutput[1]
        
        # Total Intensity
        self.I = totalI
        
