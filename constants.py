from scipy import constants as scpc
import numpy as np

mp = scpc.value(u'proton mass energy equivalent in MeV')
mn = scpc.value(u'neutron mass energy equivalent in MeV')
mN = (mp+mn)/2 # nucleon average mass
uN = (mp*mn)/(mp+mn) #nucleon reduced mass

rad2deg = 180.0/np.pi
