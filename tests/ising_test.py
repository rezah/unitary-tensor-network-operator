import pyUni10 as uni10
import itertools
import matplotlib.pyplot as plt
import matplotlib
import pylab
import numpy as np
import sys
sys.path.append('/Users/sunchong/CodingDir/TNMBL')
import qsampl 
import optuni
import random
import copy
import time
import mpo
import env
import utils
import optimize
################################## Set parameter #######################################
Model='Ising'               #could be: 'Ising' and or 'Heisenberg'
L=int(sys.argv[1])          #Numbers of particles should be even!
llayer = int(sys.argv[2])
if llayer == 1:
    L_lay = [0]
elif llayer == 2:
    L_lay = [0,1]
else:
    raise Exception("Layer number can only be 1 or 2!")
d=2                              #pysical bond-dimension  
chi=3                            #bond-dimension of MPO

J=1.0                            #coupling, 1.0
Fieldz=1.0                            #Field in z direction, 0.6


W=8                              #random interval, [-W,W]
J2=0.0                           #coupling, 0.3
hz_list=[]                        # list of randomness 
U_delta=0.0                     #if it's zero, U_list is intialize by Identity
Method='CGploy'          #methods: CGarmjo,CGploy, SVD, SteepestDescent, SteepestDescentploy 
Randomness='Fixed'
Max_number_iteratoin_SVD=60        # maximum number of sweeps for SVD method 
Max_number_iteratoin_Steepest=60   # maximum number of sweeps for SteepestDescent method
Max_number_iteratoin_CG=60   # maximum number of sweeps for SteepestDescent method
Max_SVD_iteratoin=30               #maximum number of SVD iteration
Max_Steepest_iteratoin=30          #maximum number of SteepestDescent iteration
Max_CG_iteratoin=30          #maximum number of SteepestDescent iteration
######################################################


##################MY CODE##################MY CODE##################MY CODE
##################MY CODE##################MY CODE##################MY CODE
##################MY CODE##################MY CODE##################MY CODE

try:
    sdfdsfsdf
    U_list = utils.load_ulist(L/2,llayer,outdir='./ulist/')
except:
    U_list = optuni.variantional_unitary_process(L,L_lay,Model,J,Fieldz,hz_list)
    utils.save_ulist(U_list,outdir='./ulist/')
#U_list = mpo.intialize_unitary_list( L, L_lay, d, U_delta)
#print "final: "
#print U_list

#Iteratively optimize Ulist
#U_list0 = optuni.variantional_unitary_process(L,L_lay,Model,0.1,1.,hz_list)
#U_list1 = optuni.variantional_unitary_process(L,L_lay,Model,0.5,1.,hz_list,Ulist=U_list0)
#U_list2 = optuni.variantional_unitary_process(L,L_lay,Model,J,Fieldz,hz_list,Ulist=U_list1)
    
########################################################################

def ED_spectrum():
    npId = np.eye(2)
    npSx = np.array([0.,1.,1.,0.]).reshape(2,2)
    npSy = np.array([0.,-1.,1.,0.]).reshape(2,2)
    npSz = np.array([1.,0.,0.,-1.]).reshape(2,2)
    SxSx = np.kron(npSx,npSx)

    H = np.zeros((2**L,)*2)
    for i in xrange(0,L-1):
        H += (-1.*J)*np.kron(np.kron(np.eye(2**i),SxSx),np.eye(2**(L-i-2)))
        H += (-1.*Fieldz)*np.kron(np.kron(np.eye(2**i),npSz),np.eye(2**(L-i-1)))
    H += (-1.*Fieldz)*np.kron(np.eye(2**(L-1)),npSz) # Last Z-term

    ew,ev = np.linalg.eigh(H) # non-degenerate
    np.savetxt("./data/spectrum/L%d_ed_spectrum.txt"%L, ew)
        
       

def obtain_data():

    Egs, v = qsampl.ising_ED(L,J,Fieldz)
    psi = utils.vec2uni(v,L)
    psi = qsampl.UdagPsi(U_list,psi,L)
    ave,err = qsampl.measure_pstr(U_list,L,psi,Model,J,Fieldz,coutdir='./data/coefs/',eoutdir='./data/errors/',tol=1e-1)
    print "True ED: ", Egs[0]
    print "ED from Paulis: ", ave
    print "qsample err: ", err
    
    oerr = qsampl.ising_variance(L,J,Fieldz)
    print "original error: ", oerr

#spect = qsampl.get_spectrum(U_list,Model,L,outdir="./data/spectrum/")
#ED_spectrum()
#print spect
#diag_ = UHU_spectrum()
#print np.linalg.norm(diag_-spect)

obtain_data()




#def UHU_spectrum():
#    # L=8
#    if (L != 8):
#        return
#    npId = np.eye(2)
#
#    npSx = np.array([0.,1.,1.,0.]).reshape(2,2)
#    npSy = np.array([0.,-1.,1.,0.]).reshape(2,2)
#    npSz = np.array([1.,0.,0.,-1.]).reshape(2,2)
#    SxSx = np.kron(npSx,npSx)
#
#    H = np.zeros((2**L,)*2)
#    for i in xrange(0,L-1):
#        H += (-1.*J)*np.kron(np.kron(np.eye(2**i),SxSx),np.eye(2**(L-i-2)))
#        H += (-1.*Fieldz)*np.kron(np.kron(np.eye(2**i),npSz),np.eye(2**(L-i-1)))
#    H += (-1.*Fieldz)*np.kron(np.eye(2**(L-1)),npSz) # Last Z-term 
#    
#    Hmat = utils.mat2uni(H,2**8)
#    
#    bdi = uni10.Bond(uni10.BD_IN, 2)
#    bdo = uni10.Bond(uni10.BD_OUT, 2)
#    bonds = [bdi]*8+[bdo]*8
#    Hmpo = uni10.UniTensor(bonds)
#    Hmpo.putBlock(Hmat)
#    label = []
#    for i in xrange(4):
#        s = i*2
#        W1 = 2*L
#        W2 = 4*L
#        U0 = copy.copy(U_list[i][0])
#        U0.setLabel([0,1,2,3])
#        U0T = copy.copy(U0)
#        U0T.permute([2,3,0,1],2)
#        U0T.setLabel([s+W1, s+W1+1, s,    s+1])
#        U0.setLabel( [s+L,  s+L+1,  s+W2, s+W2+1])
#        Hmpo = (U0T*Hmpo)*U0
#        label.append([s+W1, s+W1+1, s+W2, s+W2+1])
#    label = list(itertools.chain(*label))
#    label.sort()
#    Hmpo.permute(label,L)
#    Hmpo.setLabel(range(L*2))
#
#    label1 = []
#    label2 = []
#    for i in xrange(3):
#        s = i*2+1
#        U1 = copy.copy(U_list[i][1])
#        U1T = copy.copy(U1)
#        W1 = 2*L
#        W2 = 4*L
#        U1T.permute([2,3,0,1],2)
#        U1T.setLabel([s+W1, s+W1+1, s,    s+1])
#        U1.setLabel( [s+L,  s+L+1,  s+W2, s+W2+1])
#        Hmpo = (U1T*Hmpo)*U1
#        label1.append([s+W1, s+W1+1])
#        label2.append([s+W2, s+W2+1])
#    label1 = list(itertools.chain(*label1))
#    label2 = list(itertools.chain(*label2))
#    label1.sort()
#    label2.sort()
#    label = [0] + label1 + [7] + [8] + label2 + [15]
#    Hmpo.permute(label,L)
#    Hmpo.setLabel(range(L*2))
#    diag = utils.get_uni10diag(Hmpo,2**8)
#    diag = np.sort(diag)
#    print diag/8.
#    return diag
#   
# 
