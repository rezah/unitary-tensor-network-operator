import numpy as np
import sys
sys.path.append("/Users/sunchong/CodingDir/TNMBL")
from qsampl import *
from pauli2body import *
import pyUni10 as uni10
from scipy import linalg as sla
import copy 


si = np.eye(2)
sx = np.array([0.,1.,1.,0.]).reshape(2,2)
sy = np.array([0.,-1.,1.,0.]).reshape(2,2)
sz = np.array([1.,0.,0.,-1.]).reshape(2,2)

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

#def test_v210(L):
#    _,v = ising_ED(L)
#    vt = vector2uni10(v,L)
#    print v
#    print vt.getBlock()
    
        


if __name__ == "__main__":
    ind = int(sys.argv[1])
    #test_ising_ED_eigh(ind)
    test_v210(ind)


