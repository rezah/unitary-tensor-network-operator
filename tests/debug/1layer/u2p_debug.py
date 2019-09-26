import numpy as np
import sys
sys.path.append("/Users/sunchong/CodingDir/TNMBL")
from qsampl import *
from pauli2body import *
import pyUni10 as uni10
from scipy import linalg as sla
import copy 


L=4
bond_dim = 2

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

def pauli2b():
    p1b = [si,sx,sy,sz]
    p2b = []
    for i in xrange(4):
        for j in xrange(4):
            p2b.append(np.kron(p1b[i],p1b[j]))
    return p2b

bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
bdo = uni10.Bond(uni10.BD_OUT, bond_dim)
UnitU = uni10.UniTensor([bdi,bdi,bdo,bdo])
UnitU.putBlock(uni10.otimes(Si,Si))

#Ulist = [[U11,U21],[U12,U22],[U13,U23]]
Ulist = []
for i in range(2):
    Ulist.append([uni10.UniTensor('../../ulist/U%i0'%i)])

ulist = []
for i in range(2):
    ulist.append([uni2mat(Ulist[i][0], 4)])
#Ulist[1][0] = UnitU
#Ulist[0][0] = UnitU

#######################TESTS########################

def test_uhu():
    UHU = contract_UHU(L,Ulist,'Ising')
    uhu = []
    tmp = np.dot(ulist[0][0].T, np.dot(hl, ulist[0][0]))
    uhu.append(tmp)
    uhu.append(hm)
    tmp = np.dot(ulist[1][0].T, np.dot(hr, ulist[1][0]))
    uhu.append(tmp)

    for i in xrange(len(uhu)):
        print '===================='
        print 'numpy:\n', uhu[i]
        print 'mpo:\n',UHU[i].getBlock()
   
#test_uhu()
    
def test_uhu2p():
    UHU = contract_UHU(L,Ulist,'Ising')
    uhu = []
    tmp = np.dot(ulist[0][0].T, np.dot(hl, ulist[0][0]))
    uhu.append(tmp)
    uhu.append(hm)
    tmp = np.dot(ulist[1][0].T, np.dot(hr, ulist[1][0]))
    uhu.append(tmp)

    p2b = pauli2b()
    coef = []
    for l in xrange(len(uhu)):
        c_ = np.zeros(16)
        for i in xrange(16):
            c_[i] = np.trace(np.dot(p2b[i],uhu[l]))/4.
        coef.append(c_)

    Coef = UHU2pauli(L,Ulist,'Ising')

    for l in xrange(len(coef)):
        print "Numpy: \n",coef[l]
        print "MPO: \n", Coef[l]

#test_uhu2p()

def test_upsi():
    _, v = ising_ED(L)
    Psi = vec2uni10(v,L)
    Psi = UdagPsi(Ulist,Psi,L)
    v = v.reshape(2,2,2,2)
    u00 = (ulist[0][0].T).reshape(2,2,2,2)
    u10 = (ulist[1][0].T).reshape(2,2,2,2)
    v = np.einsum('ijkl,klmn->ijmn', u00,v)
    v = np.einsum('ijkl,mnkl->mnij', u10,v)
    v = v.reshape(16)
    
    print "MPO: \n"
    print Psi.getBlock()
    print "numpy: \n",
    print v

#test_upsi()

def test_psiHpsi():
    ew, v = ising_ED(L)
    E_ed = ew[0]
    psi = vec2uni10(v,L)
    psiT = copy.copy(psi)
    psiT.permute(range(L),0)
    #UHU = contract_UHU(L,Ulist,'Ising')
    Empo = 0.
    for n in range(L-1):
        H = get_Hbond_tensor('Ising',n,L)
        labelpsiT = range(4,5+n-1) + range(2)   + range(5+n-1,L+2)
        labelpsi  = range(4,5+n-1) + range(2,4) + range(5+n-1,L+2)
        psi.setLabel(labelpsi)
        psiT.setLabel(labelpsiT)
        Empo += ((psiT*H)*psi).getBlock().sum()

    print "ED: ", E_ed
    print "MPO: ", Empo


def test_psiUHUpsi():
    ew, v = ising_ED(L)   
    E_ed = ew[0]
    psi = vec2uni10(v,L)
    psi =  UdagPsi(Ulist,psi,L)
    psiT = copy.copy(psi)
    psiT.permute(range(L),0)
    UHU = contract_UHU(L,Ulist,'Ising')
    Empo = 0.
    for n in range(L-1):
        H = UHU[n] 
        labelpsiT = range(4,5+n-1) + range(2)   + range(5+n-1,L+2)
        labelpsi  = range(4,5+n-1) + range(2,4) + range(5+n-1,L+2)
        psi.setLabel(labelpsi)
        psiT.setLabel(labelpsiT)
        Empo += ((psiT*H)*psi).getBlock().sum()
    print "ED: ", E_ed
    print "MPO: ", Empo

test_psiUHUpsi()

    


def test_psiPpsi():
    pass


   





if __name__ == '__main__':
    pass
    #test_uhu2pauli()
    #test_pauli4b(ind)
    
    
    
 
