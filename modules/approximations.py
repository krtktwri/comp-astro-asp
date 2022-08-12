# curve-fit() function imported from scipy

from numpy import arccos, cos, exp

'''
Note: Input R in r_G normalized terms
'''

def beloApprox(alpha, R):
    return arccos(1-((1-cos(alpha))/(1-(1/R))))

def beloApproxCorrected(alpha, R):
    A, B, C = 0.0026091626179161505, 8.896382706025112e-06, 5.350133183173007
    return arccos(1-((1-cos(alpha))/(1-(1/R)))) - (A + B*exp(C*alpha))

