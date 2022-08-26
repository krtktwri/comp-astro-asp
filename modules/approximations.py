from numpy import arccos, cos, exp

'''
Note: Input R in r_G normalized terms
'''

def alpha2psi(alpha, R):
    return arccos(1-((1-cos(alpha))/(1-(1/R))))

def alpha2psiCorrected(alpha, R):
    A, B, C = 0.0026091626179161505, 8.896382706025112e-06, 5.350133183173007
    return arccos(1-((1-cos(alpha))/(1-(1/R)))) - (A + B*exp(C*alpha))

def psi2alpha(psi, R):
    return arccos(1 - ((1-cos(psi)) * (1-1/R)))

def psi2alphaCorrected(psi, R):
    A, B = 0.00039249005934231736, 1.6986730233029856
    return arccos(1 - ((1-cos(psi)) * (1-1/R))) + (A*exp(B*psi))
