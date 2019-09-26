import numpy as np
import sys
sys.path.append("/Users/sunchong/CodingDir/TNMBL")
from qsampl import *
from pauli2body import *
import pyUni10 as uni10
from scipy import linalg as sla
import itertools
import copy 


L=6

# ### bond construction
sx = np.array([0.,1.,1.,0.]).reshape(2,2)
sy = np.array([0.,-1,1,0.]).reshape(2,2)
sz = np.array([1.,0.,0.,-1.]).reshape(2,2)
si = np.eye(2)
Sx = uni10.Matrix(2,2,[0.,1.,1.,0.])
Sy = uni10.Matrix(2,2,[0.,-1,1,0.])
Sz = uni10.Matrix(2,2,[1.,0.,0.,-1.])
Si = uni10.Matrix(2,2,[1.,0.,0.,1.])


hm = -1.* np.kron(sx,sx)-0.5*(np.kron(sz,si)+np.kron(si,sz))
hl = -1.* np.kron(sx,sx)-0.5*np.kron(si,sz)-np.kron(sz,si)
hr = -1.* np.kron(sx,sx)-np.kron(si,sz)-0.5*np.kron(sz,si)

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


npp1b = [si,sx,sy,sz]

def np_pauli3b():
    p3b = []
    for i in xrange(4):
      for j in xrange(4):
        for k in xrange(4):
          p3b.append(np.kron(np.kron(npp1b[i],npp1b[j]),npp1b[k]))
    return p3b

def np_pauli4b():
    p4b = []
    for i in xrange(4):
      for j in xrange(4):
        for k in xrange(4):
          for l in xrange(4):
            p4b.append(np.kron(np.kron(np.kron(npp1b[i],npp1b[j]),npp1b[k]),npp1b[l]))
    return p4b

def expand_pauli(p,ind,L):
    if ind == 0:
        return np.kron(p,np.eye(2**(L-3))) 
    elif ind == (L-2):
        return np.kron(np.eye(2**(L-3)),p)
    else:
        return np.kron(np.kron(np.eye(2**(ind-1)),p),np.eye(2**(L-ind-3)))
        

pauli3b = pauli3body()
pauli4b = pauli4body()
npp3b = np_pauli3b()
npp4b = np_pauli4b()

#################################################################################

def test1():
    diag = np.zeros(2**L)
    UHU = contract_UHU(L,Ulist)   
    
    M = UHU[0].getBlock()
    l3 = 8
    l4 = 16
    m0 = np.zeros((l3,l3))
    for i in range(l3):
        for j in range(l3):
            m0[i,j] = M[i*l3+j]
    H_ = np.kron(m0,np.eye(2**(L-3)))
    diag += np.diag(H_)
    
    for n in range(1,L-2):
        M = UHU[n].getBlock()
        m0 = np.zeros((l4,l4))
        for i in range(l4):
            for j in range(l4):
                m0[i,j] = M[i*l4+j]
        H_ = np.kron( np.eye(2**(n-1)), np.kron(m0,np.eye(2**(L-3-n))))
        diag += np.diag(H_)
    
    M = UHU[-1].getBlock()
    m0 = np.zeros((l3,l3))
    for i in range(l3):
        for j in range(l3):
            m0[i,j] = M[i*l3+j]
    H_ = np.kron(np.eye(2**(L-3)),m0)
    diag += np.diag(H_)
    
    print diag


    



