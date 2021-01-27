import angmom as angmom
from itertools import groupby
from operator import itemgetter

def setup_NN_basis(j2min, j2max, tzmin, tzmax):
    basis = []
    for tz in range(tzmin,tzmax+1,1):
        for J in range(j2min,j2max+1,1):
            for S in range(0,2):
                for L in range(abs(J-S),J+S+1,1):
                    for T in range(abs(tz),2,1):
                        if ((L+S+T)%2 != 0):
                            basis_state = {}
                            basis_state['tz'] = tz
                            basis_state['l']  = L
                            basis_state['pi'] = (-1)**L
                            basis_state['s']  = S
                            basis_state['j']  = J
                            basis_state['t']  = T
                            print(basis_state)
                            basis.append(basis_state)

    print('len(basis) = ',len(basis))
    return basis

def setup_NN_channels(basis):

    states = []
    
    for bra in basis:
        for ket in basis:
            
            if angmom.kroenecker_delta(bra,ket,'j','tz','s','pi'):
                
                state = {}
                
                state['l']  = bra['l']
                state['ll'] = ket['l']
                
                state['s']  = bra['s']
                state['j']  = bra['j']
                state['t']  = bra['t']
                state['tz'] = bra['tz']
                state['pi'] = bra['pi']
                states.append(state)

    grouper = itemgetter("s", "j", "tz", "pi")
    NN_channels = []
    
    for key, grp in groupby(sorted(states, key = grouper), grouper):
        NN_channels.append(list(grp))
        

    for chn_idx, chn in enumerate(NN_channels):
        for block in chn:
            block.update({"chn_idx":chn_idx})
            
    for idx, chn in enumerate(NN_channels):
        print('channel %d: J=%d S=%d Tz=%d Pi=%d'%(idx,chn[0]['j'],chn[0]['s'],chn[0]['tz'],chn[0]['pi']))
        for block in chn:
            print('  block:',block)
    
    print('len(NN_channels) = ',len(NN_channels))
    return NN_channels
            
def lookup_NNchannel_idx(NN_channels, **kwargs):
    matching_indices = []
    for idx, chn in enumerate(NN_channels):
        for block in chn:
            if (kwargs.items() <= block.items()):
                matching_indices.append(idx)

    matching_indices = list(dict.fromkeys(matching_indices))

    return matching_indices
                
    
