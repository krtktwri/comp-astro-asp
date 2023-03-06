import numpy as np

# create a neutron star object
'''
Star class has the following attributes:
Radius (in units of Scwarzschild radius), mass (in units of solar masses), 
eta (angle between magnetic axis and spin axis), i (angle between LOS and spin axis),
B_p (magnetic field strength at the poles in 10^12 Gauss)
'''
class star:
  def __init__(self, R, M, eta, i, B_p):
    self.radius = R  # in terms of schwarzschild radius
    self.mass = M    # in terms of solar masses
    self.eta = np.deg2rad(eta)   # angle between magnetic axis and spin axis
    self.i = np.deg2rad(i)       # angle between LOS and spin axis
    self.B_p = B_p   # magnetic field strength at the poles (in 10^12 Gauss)
    
    
# Changing Angle between Magnetic Axis and Spin Axis
NS1 = star(3, 1.6, 0, 0, 1.4)
NS2 = star(3, 1.6, 45, 0, 1.4)
NS3 = star(3, 1.6, 90, 0, 1.4)

NS4 = star(3, 1.6, 0, 30, 1.4)
NS5 = star(3, 1.6, 45, 30, 1.4)
NS6 = star(3, 1.6, 90, 30, 1.4)

NS7 = star(3, 1.6, 0, 60, 1.4)
NS8 = star(3, 1.6, 45, 60, 1.4)
NS9 = star(3, 1.6, 90, 60, 1.4)

NS10 = star(3, 1.6, 0, 90, 1.4)
NS11 = star(3, 1.6, 45, 90, 1.4)
NS12 = star(3, 1.6, 90, 90, 1.4)

# Changing Observation Spectra
NS13 = star(3, 1.6, 40, 50, 1)

# Changing Magnetic Field Strength
NS14 = star(3, 1.6, 45, 60, 0.1)
NS15 = star(3, 1.6, 45, 60, 1)
NS16 = star(3, 1.6, 45, 60, 10)

# and so on