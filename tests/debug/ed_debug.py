import numpy as np
import sys
sys.path.append("/Users/sunchong/CodingDir/TNMBL")
from qsampl import *
from pauli2body import *
import pyUni10 as uni10
from scipy import linalg as sla
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
u13 = np.einsum('ij,kl->ikjl',sy,sz)
u14 = np.einsum('ij,kl->ikjl',sy,sz)
u21 = np.einsum('ij,kl->ikjl',sz,sx)
u22 = np.einsum('ij,kl->ikjl',sz,si)
u23 = np.einsum('ij,kl->ikjl',sx,sx)
npulist = [[u11,u21],[u12,u22],[u13,u23]]



U11 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U11.putBlock(uni10.otimes(Sx,Sz))
U12 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U12.putBlock(uni10.otimes(Si,Sx))
U13 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U13.putBlock(uni10.otimes(Sy,Sz))
U14 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U14.putBlock(uni10.otimes(Sy,Sz))
U21 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U21.putBlock(uni10.otimes(Sz,Sx))
U22 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U22.putBlock(uni10.otimes(Sz,Si))
U23 = uni10.UniTensor([bdi,bdi,bdo,bdo])
U23.putBlock(uni10.otimes(Sx,Sx))



Ulist = [[U11,U21],[U12,U22],[U13,U23]]

def transmat(u):
    t = u.reshape(4,4)
    return (t.T).reshape(2,2,2,2)


def test_ising_ED_getH(ind):
    '''
    Need to print H in ising_ED() function.
    '''
    if ind == 0:
        L = 2
        ising_ED(L)
        H = -1.*np.kron(sx,sx)
        H += -1.*np.kron(si,sz)
        H += -1.*np.kron(sz,si)
        
        print "Hand made:"
        print H
    if ind == 1:
        L = 4
        ising_ED(L)
        H  = -1.*np.kron(np.kron(sx,sx),np.kron(si,si))
        H += -1.*np.kron(np.kron(si,np.kron(sx,sx)),si)
        H += -1.*np.kron(np.kron(si,si),np.kron(sx,sx))
        H += -1.*np.kron(np.kron(sz,np.kron(si,si)),si)
        H += -1.*np.kron(np.kron(si,np.kron(sz,si)),si)
        H += -1.*np.kron(np.kron(np.kron(si,si),sz),si)
        H += -1.*np.kron(np.kron(si,si),np.kron(si,sz))
        print "Hand Made: "
        print H
    if ind == 2:
        L = 14 # raise error
        ising_ED(L)
    if ind ==3:
        L = 3 # raise error
        ising_ED(L) 

def test_ising_ED_eigh(L):
    ising_ED(L) 

def test_udagpsi():

    print npulist[2][0].reshape(4,4)
    E, v = ising_ED(L)
    psi = vector2uni10(v,L)
    UdagPsi(Ulist,psi,L)
    # First layer
    v = v.reshape([2,]*L)
    v = np.einsum('ijklmn,pqij->pqklmn',v,transmat(npulist[0][0]))
    v = np.einsum('klijmn,pqij->klpqmn',v,transmat(npulist[1][0]))
    v = np.einsum('klmnij,pqij->klmnpq',v,transmat(npulist[2][0]))

    # second layer
    v = np.einsum('kijlmn,pqij->kpqlmn',v,transmat(npulist[0][1]))
    v = np.einsum('klmijn,pqij->klmpqn',v,transmat(npulist[1][1]))
    v = v.reshape(2**L)
    
    print v
    #TODO test with numpy
    #i = 0
    #U14_ = copy.copy(U14)
    #U14_.permute([2,3,0,1],2)
    #U14_.setLabel([i,i+1,i+L,i+L+1])
    #print U14
    #print U14_
    
   
        


if __name__ == "__main__":
    ind = int(sys.argv[1])
    #test_ising_ED_eigh(ind)
    test_udagpsi()


