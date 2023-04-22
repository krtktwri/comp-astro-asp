from numpy import arccos, cos, exp, deg2rad, rad2deg, arange
from modules.converter import alpha2l, theta
from modules.exact import trajectory
from modules.myconstants import r_G

'''
Note: Input R in r_G normalized terms
'''

def alpha2psi(alpha, R):
    return arccos(1-((1-cos(alpha))/(1-(1/R))))

def alpha2psiCorrected(alpha, R):
    A, B, C = 0.0026091626179161505, 8.896382706025112e-06, 5.350133183173007
    return arccos(1-((1-cos(alpha))/(1-(1/R)))) - (A + B*exp(C*alpha))

def psi2alpha(psi, R):
    return arccos(1/R + cos(psi)*(1-1/R))

def psi2alphaCorrected(psi, R):
    A, B = 0.00039249005934231736, 1.6986730233029856
    return arccos(1 - ((1-cos(psi)) * (1-1/R))) + (A*exp(B*psi))

'''
Function that ouputs exact relation for a given R between alpha and phi
'''
def alpha_for_R(R):
    
    # Initializing lists and Arrays
    alphaArray = arange(deg2rad(5), deg2rad(89.5), deg2rad(1))
    sols, lList = [], []

    # Converting alpha values to impact parameters
    for a in alphaArray:
        lList.append(alpha2l(a, R))

    # Calculating trajectories
    for i in range(0, len(lList)):
        sols.append(trajectory(lList[i], R*r_G, False))

    # extracting the angle at which the ray emerges    
    psi = [] # Angle at which the ray emerges
    for i in range(0, len(sols)):
        psi.append(theta(sols[i][0][-1], sols[i][1][-1]))

    return alphaArray, psi

'''
Function to compute a fourth order polynomial given the independent variable and the polynomial coefficients
This is used to fit the residues of beloborodov as well as to fit the trend of the coefficients used to fit through beloborodov
'''
def fourth_order_poly(x, A, B, C, D):             
    return  A*x**3 + B*x**2 + C*x + D

'''
The modified beloborodov approximation being employed in computations  
'''

# Coefficients for a paramteric fit through correction coefficients
A1, A2, A3, A4 = [-1.382364238031976228e-08,1.261611907346646713e-07,-4.007535578418975164e-07,4.589947962920749492e-07]
B1, B2, B3, B4 = [3.179392984048052565e-06,-2.823471312371944979e-05,8.603088610468706348e-05,-9.240985994226368126e-05]
C1, C2, C3, C4 = [-1.926743652254383446e-04,1.681248251252310394e-03,-4.976838196856497171e-03,5.136783637201085967e-03]
D1, D2, D3, D4 = [2.833930597530380773e-03,-2.439343625511589814e-02,7.110382899420360225e-02,-7.019667298595769211e-02]

'''
Following function takes Radius as input and outputs the coefficients 
of the fourth order polynomial that fits through the residues of Beloborodov approximation
'''
def extract_coefficients(R):
    
    a = fourth_order_poly(R, A1, A2, A3, A4)
    b = fourth_order_poly(R, B1, B2, B3, B4)
    c = fourth_order_poly(R, C1, C2, C3, C4)
    d = fourth_order_poly(R, D1, D2, D3, D4)
    
    return a, b, c, d

'''
Following function computes, to a greater accuracy than belo2003, the emission angle Alpha
from the input phi (i. e. latitude of the emission point), radius (expressed in units of gravitational radius) 
and the correction coefficients (corr_coeff = extract_coefficients(R))
'''

def fastApprox(phi, R, corr_coeff):
    
    a, b, c, d = corr_coeff
    residue = fourth_order_poly(rad2deg(phi), a, b, c, d)
    alpha = psi2alpha(phi, R) + residue
    
    return alpha, residue