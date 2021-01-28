from NN_potential import chiral_LO as potential

# for importing Cython chiral potential. DO NOT USE
#import sys
#sys.path.append('/Users/andreas/software/chp/chp_python')
#import pychp as pychp
import angmom as angmom
import numpy as np
import constants as const

# <pL|V|ppLL>_SJTtz
def V(p,pp,L,LL,S,J,T,Tz):

    # fac: pi factor (and sign-convention in L-LL-off-diagonal coupled channels)
    # RM corresponds to fac = +1 for idx 4,5 
    
    assert(angmom.triag(L,S,J))
    assert(angmom.triag(LL,S,J))
    if L == LL:
        fac = +np.pi/2.0
        if L<J:
            # --
            coup = 1
            idx  = 3
        elif L>J:
            # ++
            coup = 1
            idx  = 2
        else:
            if S==1:
                coup = 0
                idx  = 1
            else:
                coup = 0
                idx  = 0
    else:
        if L<J:
            # -+
            coup = 1
            idx  = 5
            fac  = -np.pi/2
        else:
            # +-
            coup = 1
            idx  = 4
            fac  = -np.pi/2

    # Seans potential
    # [V_uncoupled_S0, V_uncoupled_S1, V_coupled_pp, V_coupled_mm, V_coupled_pm, V_coupled_mp]
    V = potential.V(p,pp,coup,S,J,T,Tz)[idx]*fac

    # pychp
    #V = pychp.V(p,pp,coup,S,J,T,Tz)[idx]*fac
    
    return V

def Vmtx(pmesh,L,LL,S,J,T,Tz):
    mtx = np.zeros((len(pmesh), len(pmesh)))
    for pidx, p in enumerate(pmesh):
        for ppidx, pp in enumerate(pmesh):
            mtx[pidx][ppidx] = V(p,pp,L,LL,S,J,T,Tz)
    return np.array(mtx)

def Vvec(pmesh,po,L,LL,S,J,T,Tz,direction='col'):
    vec = np.zeros(len(pmesh))
    if direction is 'col':
        for pidx, p in enumerate(pmesh):
            vec[pidx] = V(p,po,L,LL,S,J,T,Tz)
    elif direction is 'row':
        for pidx, p in enumerate(pmesh):
            vec[pidx] = V(po,p,L,LL,S,J,T,Tz)

    return np.array(vec)

def init_pychp():

    potential_model = pychp.preset_lo
    #potential_model = pychp.preset_idaho_n3lo
    #potential_model = pychp.preset_nnlosat
    #potential_model = pychp.preset_nnloopt
    #potential_model = pychp.preset_dnnlo
    
    # LO lecs for Lambda 450
    lambda_cutoff = 450.0 # MeV
    lecs = np.array([[-0.112927,-0.087340]],dtype=np.double)
    
    # DNNLO_go(394) lecs
    #lambda_cutoff = 394.0 # MeV
    #lecs = np.array([[-0.338142,-0.339250,-0.338746,-0.259839,+2.5053890,
    #     +0.7004990,-0.3879600,-0.9648560,+1.0021890,+0.4525230,
    #     -0.8883122,-0.74,-0.49,-0.65,+0.96]],dtype=np.double)
    # initialize potential
    potential_model(lecs[0], lambda_cutoff)


# setup potential matrices Vchn for all channels:
# - coupled   channels: dim(Vchn) = [2Np] x [2Np]
# - uncoupled channels: dim(Vchn) = [Np]  x [Np]
# where Np = len(x) momentum mesh, in MeV/c
# split mode [True/False] determines if the potential matrices should be split along the LEC-value-axis
def setup_V(NN_channel, pmesh, Tlab):

    #init_pychp()
    
    Np = len(pmesh)
    
    nof_chns = len(NN_channel)
    Np_chn = np.sqrt(nof_chns)*Np
        
    chn_idx = NN_channel[0]['chn_idx']
    print('populating Vmtx in channel %d'%(chn_idx))
        
    m = []

    for idx, block in enumerate(NN_channel):
            
        l  = block['l']
        ll = block['ll']
        s  = block['s']
        j  = block['j']
        t  = block['t']
        tz = block['tz']
        
        if NN_channel[0]['tz'] == -1:
            mu = const.mp/2
            ko2 = 2*const.mp*Tlab
        elif NN_channel[0]['tz'] ==  0:
            mu = const.uN
            ko2 = const.mn**2*Tlab*(Tlab+2*const.mp)/((const.mp+const.mn)**2 + 2*Tlab*const.mn)
        elif NN_channel[0]['tz'] == +1:
            mu = const.mN/2
            ko2 = 2*const.mp*Tlab
        else:
            print('unknown isospin')
            
        ko = np.sqrt(ko2)
                
        pmesh1 = np.hstack((pmesh,[ko]))
        mtx = np.copy(Vmtx(pmesh1,l,ll,s,j,t,tz))

        m.append(mtx)

    if len(NN_channel) >1:
        V = np.copy(np.vstack((np.hstack((m[0],m[1])),
                                  np.hstack((m[2],m[3])))))
    else:
        V = np.copy(m[0])
                   
    return V, ko
