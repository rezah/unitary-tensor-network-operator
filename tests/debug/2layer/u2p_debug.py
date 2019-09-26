import numpy as np
import sys
sys.path.append("/Users/sunchong/CodingDir/TNMBL")
from qsampl import *
#from pauli2body import *
import pyUni10 as uni10
from scipy import linalg as sla
import copy 
import utils


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

#test_ulist2p()
#test_u2p()
    
def test_pauli_str_contract():
    strA = '2123'
    strB = '1012'
    posA = [0,1]
    posB = [2,3]
    nstr,coef = pauli_str_contract(strA,strB,posA,posB)
    print 'correct str: ', '103323', ', correct coef:', '-1'
    print nstr, coef

    
def test_pauli3b(ind):
    p3b = pauli3body()
    if ind == 0:
        print '------------001------------'
        print np.kron(np.kron(si,si),sx)
        print p3b[1].getBlock()
    if ind ==1:
        print '------------012------------'
        print np.kron(np.kron(si,sx),sy)
        print p3b[6].getBlock()
    if ind ==2:
        print '------------123------------'
        print np.kron(np.kron(sx,sy),sz)
        print p3b[27].getBlock()

def test_pauli4b(ind):
    p4b = pauli4body()
    if ind == 0:
        print '-----------0001------------'
        print np.kron(np.kron(np.kron(si,si),si),sx)
        print p4b[1].getBlock()
    if ind == 1:
        print '-----------1021------------'
        print np.kron(np.kron(np.kron(sx,si),sy),sx)
        print p4b[73].getBlock()
    if ind == 2:
        print '-----------0123------------'
        print np.kron(np.kron(np.kron(si,sx),sy),sz)
        print p4b[27].getBlock()



def pauli3b():
    p1b = [si,sx,sy,sz]
    p3b = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                p3b.append(np.kron(np.kron(p1b[i],p1b[j]),p1b[k]))

    return p3b

def pauli4b():
    p1b = [si,sx,sy,sz]
    p4b = []
    for i in range(4):
      for j in range(4):
        for k in range(4):
          for l in range(4):
            p4b.append(np.kron(np.kron(np.kron(p1b[i],p1b[j]),p1b[k]),p1b[l]))

    return p4b

    
def test_uhu2pauli(ind):
    pcoef = UHU2pauli(L,Ulist)
    p3b = pauli3b()
    p4b = pauli4b()
    if ind == 0:
        u11_ = u11.reshape(4,4)
        uhu1 = np.dot(u11_.T,np.dot(hl,u11_))
        uhu1 = uhu1.reshape(2,2,2,2)
        u21T = (u21.reshape(4,4).T).reshape(2,2,2,2)
        uhu = np.einsum('ijkl,mknp,plrs->mijnrs',u21T,uhu1,u21).reshape(8,8)
        npcoef = []
        for i in range(64):
            cnp = np.sum(np.diag(np.dot(uhu,p3b[i])))/8.
            npcoef.append(cnp)
        print pcoef[0]
        print npcoef

    if ind == 1:
        hm_ = hm.reshape(2,2,2,2)
        u11T = (u11.reshape(4,4).T).reshape(2,2,2,2)
        u12T = (u12.reshape(4,4).T).reshape(2,2,2,2)
        u21T = (u21.reshape(4,4).T).reshape(2,2,2,2)
        uhu1 = np.einsum('ijkl,mnpq,lprs,kruv,sqwx->ijmnuvwx',u11T,u12T,hm_,u11,u12)
        uhu  = np.einsum('ijkl,mklnpqrs,qruv->mijnpuvs',u21T,uhu1,u21).reshape(16,16)
        npcoef = []
        for i in range(256):
            cnp = np.trace(np.dot(uhu,p4b[i]))/16.
            npcoef.append(cnp)
        print pcoef[1]
        print npcoef

    if ind == 2:
        u12_ = u12.reshape(4,4)
        uhu1 = np.dot(u12_.T,np.dot(hm, u12_)).reshape(2,2,2,2)
        u21T = (u21.reshape(4,4).T).reshape(2,2,2,2)
        u22T = (u22.reshape(4,4).T).reshape(2,2,2,2)
        uhu = np.einsum('ijkl,mnpq,lpxy,kxrs,yquv->ijmnrsuv',u21T,u22T,uhu1,u21,u22).reshape(16,16)        
        npcoef = []
        for i in range(256):
            cnp = np.trace(np.dot(uhu,p4b[i]))/16.
            npcoef.append(cnp)
        print pcoef[2]
        print npcoef

    if ind == 3:
        u13_ = u13.reshape(4,4)
        uhu1 = np.dot(u13_.T,np.dot(hr,u13_)).reshape(2,2,2,2)
        u22T = (u22.reshape(4,4).T).reshape(2,2,2,2)
        uhu  = np.einsum('ijkl,lmqn,kqrs->ijmrsn',u22T,uhu1,u22).reshape(8,8)        
        npcoef = []
        for i in range(64):
            cnp = np.sum(np.diag(np.dot(uhu,p3b[i])))/8.
            npcoef.append(cnp)
        print pcoef[-1]
        print npcoef
        


        






if __name__ == '__main__':
    ind = int(sys.argv[1])
    #test_uhu2pauli()
    #test_pauli4b(ind)
    test_uhu2pauli(ind)
    
    
    
 
