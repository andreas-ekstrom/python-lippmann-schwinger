import numpy as np

class OPEP():
	
	def __init__(self):
		
		# Create Gauss-Legendre grid to perform angular integration
		angLattice = np.polynomial.legendre.leggauss(96)
		self.x = angLattice[0]			  
		self.w = angLattice[1]
		
		self.gA 		= 1.289 	 # axial coupling constant 			  [no units]
		self.fpi 		= 92.2		 # pion decay constant 				  [MeV]        (used convention 130.41/sqrt(2) = 92.2)
		self.mpi 		= 138.039 	 # pion_0 mass 					  [MeV]        (average of pi+, pi-, and pi0 masses)

	def find_root(self, d, x):
		#using recurrance relation for Legendre polynomials
		if d==0:
			return 1.
		elif d==1:
			return x
		else:
			Pnext = 0
			Pprev = 1
			Pcurr = x
			c = 1
			while c<d:
				Pnext = ((2*c+1)*x*Pcurr - c*Pprev)/float(c+1)	#recurrence relation
				Pprev = Pcurr
				Pcurr = Pnext
				c += 1
			return Pnext

	def isoFac(self, L, S):
		
		# if L+S is even then T must be odd (1), and vice versa
		
		T = (1-L-S & 1) #0 if even, 1 if odd
		
		isoFactor = -3*(1-T) + 1*T	# T is either 0 or 1
		return isoFactor

	# All taken from Erkelenz's paper on "momentum space calculations and helicity formalism in nuclear physics" (Nucl. Phys. A176 (1971) 413-432)
	def V(self, qi, qo, coupled, J):
		
		# Angular integrals
		# Can improve efficienty here by using previous J calculations of integral_
		integral_0 = self.ang_Integral(qi,qo,J,0)
		integral_P = self.ang_Integral(qi,qo,J+1,0)
		integral_M = self.ang_Integral(qi,qo,J-1,0)
		
		# declare return values so they're all defined
		V_uncoupled_S0 = 0
		V_uncoupled_S1 = 0
		V_coupled_mm   = 0
		V_coupled_pm   = 0
		V_coupled_mp   = 0
		V_coupled_pp   = 0
		
		if coupled==False:
			integral_1 = self.ang_Integral(qi,qo,J,1)

			# OPEP uncoupled interactions
			V_uncoupled_S0 = 2 * (-(qo**2+qi**2)*integral_0 + 2*qo*qi*integral_1)
			V_uncoupled_S1 = 2 * ((qo**2+qi**2)*integral_0 - 2*qo*qi*(1./(2*J+1))*(J*integral_P + (J+1)*integral_M))

			V_uncoupled_S0 *= self.isoFac(J,0)
			V_uncoupled_S1 *= self.isoFac(J,1)
		else:
			# OPEP coupled interactions ( "m": l-1=J, "p": l+1=J)

			if J!=0:	# These cannot exist for 3P0-wave
				V_coupled_mm    = (2./(2*J+1)) * ((qo**2+qi**2)*integral_M - 2*qo*qi*integral_0)
				V_coupled_pm    = (4*np.sqrt(J*(J+1))/(2*J+1)) * (qi**2*integral_P + qo**2*integral_M - 2*qo*qi*integral_0)
				V_coupled_mp    = (4*np.sqrt(J*(J+1))/(2*J+1)) * (qi**2*integral_M + qo**2*integral_P - 2*qo*qi*integral_0)
			
				V_coupled_mm   *= self.isoFac(J-1,1)
				V_coupled_pm   *= self.isoFac(J-1,1)
				V_coupled_mp   *= self.isoFac(J+1,1)

			V_coupled_pp    = (2./(2*J+1)) * (-(qo**2+qi**2)*integral_P + 2*qo*qi*integral_0)
			V_coupled_pp   *= self.isoFac(J+1,1)
			
		return [V_uncoupled_S0, V_uncoupled_S1, V_coupled_mm, V_coupled_pm, V_coupled_mp, V_coupled_pp]
			
	def ang_Integral(self, qi, qo, J, l):
		#this is roughly 20x faster than writing it explicitly in a for-loop
		integrand = np.sum(self.w * self.pot_OPEP_mom(qi,qo,self.x) * self.x**l * self.find_root(J, self.x))
		
		return np.pi*integrand
	
	def pot_OPEP_mom(self, qi, qo, z):
		
		q2 = qi**2 + qo**2 - 2*qi*qo*z
	
		return -(self.gA**2/(4*self.fpi**2))*(1./(q2+self.mpi**2))
		
