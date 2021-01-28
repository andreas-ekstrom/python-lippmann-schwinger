import quantum_states as qs
import mesh as mesh
import scattering as sc
import potential as pot
import numpy as np
import auxiliary as aux

# NN basis states
basis = qs.setup_NN_basis(j2min=0,j2max=2,tzmin=0,tzmax=0)
# channels: all 2b channels with conserved J T Pi (and S)
NN_channels = qs.setup_NN_channels(basis)

# setup a gauss legendre momentum mesh in MeV/c
Np = 100

p, w = mesh.gauss_legendre_inf_mesh(Np)

Tlab = 10.0

#select a specific channel
NN_channel = NN_channels[5]

#setup the strong interaction potential in this channel
#(convert Tlab to prel)
V, ko = pot.setup_V(NN_channel,p,Tlab)

#solve the T matrix in a specific channel
T = sc.compute_Tmatrix(NN_channel,V,ko,p,w)
phase = sc.compute_phase_shifts(NN_channel,ko,T)
print(phase)
