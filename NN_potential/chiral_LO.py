
import numpy as np

from . import OPEP

Lambda	= 450	 		 # cut-off for renormalization of LO  [MeV]
C1S0	= -0.112927/100. # contact term C1S0 for lambda = 450 [MeV]
C3S1	= -0.087340/100. # contact term C3S1 for lambda = 450 [MeV]
                                         
Mp 	= 938.272 		     # proton mass  					  [MeV]
Mn 	= 939.565 		     # neutron mass 					  [MeV]
	
# V(...) returns a list of length 6
# The elements are on the form [V_S0, V_S1, V_pp, V_mm, V_pm, V_mp]
# where S0-> S=0, S1-> S=1, mm-> l=l'=J-1, mp-> l=J-1, l'=J+1, etc
def V(qi, qo, coupled, S, J, T, Tz):

	# model.V returns a list of length 6
	# The elements are on the form [V_S0, V_S1, V_pp, V_mm, V_pm, V_mp]
	# where S0-> S=0, S1-> S=1, mm-> l=l'=J-1, mp-> l=J-1, l'=J+1, etc
	V_elements = OPEP.V(qi, qo, coupled, J)
	
	# LO contact terms
	if J==0: V_elements[0] += C1S0 # contact term C1S0 at leading order
	if J==1: V_elements[3] += C3S1 # contact term C3S1 at leading order

	# define neucleon mass (as used in Machleidt)
	if Tz==0:
		MN = 2*Mn * Mp / (Mn + Mp)
	elif Tz==-1:	# use idiot-convention Tz<0 for positively-charged particle
		MN = Mp
	else:
		MN = Mn

	Epi = np.sqrt(MN**2 + qi**2) 		# relativistic energy of in-going state
	Epo = np.sqrt(MN**2 + qo**2) 		# relativistic energy of out-going state
			
	relFactor_i = np.sqrt(MN/Epi)		# minimal-relativity factor of in-going state
	relFactor_o = np.sqrt(MN/Epo)		# minimal-relativity factor of out-going state
	
	const = relFactor_o*relFactor_i

	# regulator-functions and Fourier transform constants
	f1 = np.exp(-(qi/Lambda)**6)
	f2 = np.exp(-(qo/Lambda)**6)

	const *= f1*f2/(2*np.pi)**3

	# Python list-comprehension: multiplies all elements in "V_elements" by "const"
	V_elements = [element*const for element in V_elements]

	return V_elements
