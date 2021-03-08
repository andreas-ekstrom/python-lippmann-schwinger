import quantum_states as qs
import mesh as mesh
import scattering as sc
import potential as pot
import numpy as nppyt
import auxiliary as aux
import csv

# NN basis states
basis = qs.setup_NN_basis(j2min=0,j2max=2,tzmin=0,tzmax=0)


# channels: all 2b channels with conserved J T Pi (and S)
"""
2b = two body?
https://www.asc.ohio-state.edu/physics/ntg/8805/notes/section_5_Scattering_2.pdf
"""
NN_channels = qs.setup_NN_channels(basis)

# setup a gauss legendre momentum mesh in MeV/c
#Np = 100
"""
https://numpy.org/doc/stable/reference/generated/numpy.polynomial.legendre.leggauss.html
"""


#print(p[99])
Tlab = 100.0
#select a specific channel
NN_channel = NN_channels[0]

#temp = pot.V(5,5,0,0,0,0,1,0)
#print(temp)

#setup the strong interaction potential in this channel
#(convert Tlab to prel)
with open("fasN_py.csv", "w", newline="") as file:
    writer = csv.writer(file)
    for Np in range(3,501):
        p, w = mesh.gauss_legendre_inf_mesh(Np)
        V, ko = pot.setup_V(NN_channel,p,Tlab)
    #print(V[10,10])

    #solve the T matrix in a specific channel
        T = sc.compute_Tmatrix(NN_channel,V,ko,p,w)
        phase = sc.compute_phase_shifts(NN_channel,ko,T)
        #phase = phase.real
        phase = phase[0][0]
        phase = phase.real
        print(phase)
        print(Np)
        writer.writerow([phase, Np])
        
print("klar :)")

#print(phase)
