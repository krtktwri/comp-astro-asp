from numpy import arccos, cos

'''
Note: Input R in r_G normalized terms
'''

def beloApprox(alpha, R):
    return arccos(1-((1-cos(alpha))/(1-(1/R))))