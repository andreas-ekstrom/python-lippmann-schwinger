import numpy as np
import potential as pot
import constants as const
import auxiliary as aux
from scipy import linalg

def setup_G0_vector(p, w, ko):
    
    N = len(p)

    D = np.zeros((2*N+2), dtype=complex)
    
    D[0:N] = w*p**2/(ko**2 - p**2)		# Gaussian integral

    #print('   G0 pole subtraction')
    D[N]  = -np.sum( w/(ko**2 - p**2 ) )*ko**2 	# Principal value
    D[N] -= 1j*ko * (np.pi/2)
    
    D[N+1:2*N+2] = D[0:N+1]
            
    return D

def setup_VG_kernel(chn,V,ko,p,w):

    if chn[0]['tz'] == -1:
        mu = const.mp/2
    if chn[0]['tz'] ==  0:
        mu = const.uN
    if chn[0]['tz'] == +1:
        mu = const.mN/2
            
    Np = len(p)

    nof_blocks = len(chn)
    chn_idx = chn[0]['chn_idx']
    print('setting up G_0(ko) in channel %d'%(chn_idx))
    
    Np_chn = np.int(np.sqrt(nof_blocks)*(Np+1))

    # Go-vector dim(u) = 2*len(p)+2
    G0 = setup_G0_vector(p, w, ko)

    g = np.copy(G0[0:Np_chn])
    VG = np.zeros((len(g),len(g)),dtype=complex)
            
    for g_idx, g_elem in enumerate(g):
        VG[:,g_idx] = V[:,g_idx] * g_elem * (2*mu)
                
    return VG

def compute_Tmatrix(NN_channel,V,ko,p,w):

    chn_idx = NN_channel[0]['chn_idx']
 
    VG = setup_VG_kernel(NN_channel,V,ko,p,w)
    # T = (1-VG)^{-1}V
    eye_VG =  np.eye(VG.shape[0]) - (2.0/np.pi)*VG
    print('solving for the complex T-matrix in channel %d'%(chn_idx))
    Tmtx = np.matmul(np.linalg.inv(eye_VG),V)

    return Tmtx

def blattToStapp(delta_minus_BB, delta_plus_BB, twoEpsilonJ_BB):
    
    # Stapp convention (bar-phase shifts) in terms of Blatt-Biedenharn convention
    twoEpsilonJ = np.arcsin(np.sin(twoEpsilonJ_BB)*np.sin(delta_minus_BB - delta_plus_BB))	# mixing parameter
    delta_minus	= 0.5*(delta_plus_BB + delta_minus_BB + np.arcsin(np.tan(twoEpsilonJ)/np.tan(twoEpsilonJ_BB)))*const.rad2deg
    delta_plus	= 0.5*(delta_plus_BB + delta_minus_BB - np.arcsin(np.tan(twoEpsilonJ)/np.tan(twoEpsilonJ_BB)))*const.rad2deg
    epsilon	= 0.5*twoEpsilonJ*const.rad2deg
    
    return delta_minus, delta_plus, epsilon

def compute_phase_shifts(NN_channel,ko,T):

    phases = []
    
    nof_blocks = len(NN_channel)
    chn_idx = NN_channel[0]['chn_idx']
    print('computing phase shifts in channel %d'%(chn_idx))
        
    Np = T.shape[0]
    if NN_channel[0]['tz'] == -1:
        mu = const.mp/2
    if NN_channel[0]['tz'] ==  0:
        mu = const.uN
    if NN_channel[0]['tz'] == +1:
        mu = const.mN/2
        
    # S/T matrix relation
    #S = 1-4i*mu*ko*T
            
    fac = 2*mu*ko
        
    if nof_blocks > 1:
        Np = np.int((Np-2)/2)
        T11 = T[Np,Np]
        T12 = T[2*Np+1,Np]
        T22 = T[2*Np+1,2*Np+1]
            
        # Blatt-Biedenharn (BB) convention
        twoEpsilonJ_BB = np.arctan(2*T12/(T11-T22))	# mixing parameter
        delta_plus_BB  = -0.5*1j*np.log(1 - 1j*fac*(T11+T22) + 1j*fac*(2*T12)/np.sin(twoEpsilonJ_BB))
        delta_minus_BB = -0.5*1j*np.log(1 - 1j*fac*(T11+T22) - 1j*fac*(2*T12)/np.sin(twoEpsilonJ_BB))
        
        delta_minus, delta_plus, epsilon = blattToStapp(delta_minus_BB, delta_plus_BB, twoEpsilonJ_BB)
        phases.append([delta_minus,delta_plus, epsilon])
    else:
        # uncoupled
        Np = Np-1
        T = T[Np,Np]
        Z = 1-fac*2j*T
        # S=exp(2i*delta)
        delta = (-0.5*1j)*np.log(Z)*const.rad2deg
        
        phases.append([delta])
        
    return phases
