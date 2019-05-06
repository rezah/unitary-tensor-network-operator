import numpy as np
import sys
sys.path.append("/Users/sunchong/CodingDir/TNMBL")
from qsampl import *
import pyUni10 as uni10
from scipy import linalg as sla
import copy 


L=6

# ### bond construction
sx = np.array([0.,1.,1.,0.]).reshape(2,2)
sy = np.array([0.,-1.j,1.j,0.]).reshape(2,2)
sz = np.array([1.,0.,0.,-1.]).reshape(2,2)
si = np.eye(2)
Sx = uni10.Matrix(2,2,[0.,1.,1.,0.])
Sy = uni10.Matrix(2,2,[0.,-1.j,1.j,0.])
Sz = uni10.Matrix(2,2,[1.,0.,0.,-1.])
Si = uni10.Matrix(2,2,[1.,0.,0.,1.])

bdi = uni10.Bond(uni10.BD_IN,2)
bdo = uni10.Bond(uni10.BD_OUT,2)

u11 = np.einsum('ij,kl->ikjl',sx,sz)
u12 = np.einsum('ij,kl->ikjl',si,sx)
u13 = np.einsum('ij,kl->ikjl',sz,sz)
u21 = np.einsum('ij,kl->ikjl',sz,sx)
u22 = np.einsum('ij,kl->ikjl',sz,si)
u23 = np.einsum('ij,kl->ikjl',sx,sx)


U11 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U11.putBlock(uni10.otimes(Sx,Sz))
U12 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U12.putBlock(uni10.otimes(Si,Sx))
U13 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U13.putBlock(uni10.otimes(Sz,Sz))
U21 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U21.putBlock(uni10.otimes(Sz,Sx))
U22 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U22.putBlock(uni10.otimes(Sz,Si))
U23 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U23.putBlock(uni10.otimes(Sx,Sx))

Ulist = [[U11,U21],[U12,U22],[U13,U23]]

def test_u2p():
    umat = uni10.Matrix(4,4,[0.,0.,0.9,0.,0.,0.,0.,0.1,0.9,0.,0.,0.,0.,0.1,0.,0.])
    umat = uni10.Matrix(4,4,[-1., -0., -0., -1.,-0., -0., -1., -0., -0., -1., -0., -0.,-1.,-0.,-0.,1.])
    U = uni10.UniTensor([bdi,bdi,bdo,bdo],'U')
    U.putBlock(umat)
    pstr, _ = pauli2body()
    coef = []
    for i in range(len(pstr)):
        c = unitary2pauli(U,pstr[i],2)
        coef.append(c)

    print coef

def test_ulist2p():
    Ucoef = ulist2pauli(Ulist,2)
    for i in range(len(Ucoef)):
        print Ucoef[i]

test_ulist2p()
#test_u2p()
    
    
 
