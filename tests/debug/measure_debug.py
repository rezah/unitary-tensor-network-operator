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

def test_str2int():
    mystr = [1,3,2]
    base = 4
    print str2int(mystr, base)

def test_coefanalys():
    import pickle
    pcoefs  = UHU2pauli(L,Ulist,pauli3b,pauli4b)
    print len(pcoefs[-1])
    print pcoefs[-1]
    coef_analys(pcoefs, w2file=True,datadir='./output/')
    with open("./zcoef.txt", "rb") as fp:
        b = pickle.load(fp)
    with open("./nzcoef.txt", "rb") as fp:
        c = pickle.load(fp)
    print b[-1]
    print c[-1]

def test_psipaulipsi(ind):
    Egs, v = ising_ED(L)
    coefs = UHU2pauli(L,Ulist,pauli3b,pauli4b)
    psi = vec2uni10(v,L)
    l3b = 64
    l4b = 256
    npe = []
    npe3 = np.zeros(l3b)
    for i in xrange(l3b):
        if abs(coefs[0][i]) < 1e-10:
            continue
        pall = expand_pauli(npp3b[i],0,L)
        npe3[i] = np.dot(v.T,np.dot(pall,v))
    npe.append(npe3)
    for n in range(1,L-2):
        coef4 = coefs[n]
        npe4 = np.zeros(l4b)
        for i in xrange(l4b):
            if abs(coef4[i]) < 1e-10:
                continue
            pall = expand_pauli(npp4b[i],n,L)
            npe4[i] = np.dot(v.T,np.dot(pall,v))
        npe.append(npe4)

    npe3 = np.zeros(l3b)
    for i in xrange(l3b):
        if abs(coefs[-1][i]) < 1e-10:
            continue
        pall = expand_pauli(npp3b[i],L-2,L)
        npe3[i] = np.dot(v.T,np.dot(pall,v))
    npe.append(npe3)

    ave = measure_pstr(Ulist, L, psi)

    print "numpy: "
    print npe[ind]
    print "mycode: "
    print ave[ind]

def test_measure():
    Egs, v = ising_ED(L)
    psi = vec2uni10(v,L)
    ave,err = measure_pstr(Ulist,L,psi)
    print "average: ", ave
    print "error: ", err



    
    

    


    


if __name__ == "__main__":
    ind = int(sys.argv[1])
    test_measure()
    

    
