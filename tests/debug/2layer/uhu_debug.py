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

def test_hbond(ind):

    Hm = get_Hbond_tensor(edge=False)
    Hl = get_Hbond_tensor(edge=True, edgeside='L')
    Hr = get_Hbond_tensor(edge=True, edgeside='R')
    
    
    hm = -1.* np.kron(sx,sx)-0.5*(np.kron(sz,si)+np.kron(si,sz))
    hl = -1.* np.kron(sx,sx)-0.5*np.kron(si,sz)-np.kron(sz,si)
    hr = -1.* np.kron(sx,sx)-np.kron(si,sz)-0.5*np.kron(sz,si)

    nr = [hm,hl,hr]
    mr = [Hm,Hl,Hr]
    label = ["interior", "left edge", "right edge"]
    
    print "Bond ID: ", label[ind] 
    print "Numpy result: "
    print nr[ind]
    print "My result: "
    print mr[ind]


# ### U-mpo
def test_Umpo(ind,dag=False):

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
    
    if ind == 0:
        print "Bond ID: Left Bond"
    
        U = get_unitary4bond_up(Ulist,L,0)
        u = np.einsum('ijkl,lmnp->ijmknp',u11,u21).reshape(8,8)
        if dag is False:
            print "Numpy result: "
            print u
            print "My result: "
            print U.getBlock()
        else:
            U.permute([5,6,7,2,3,4],3)
            U.setLabel([8,9,10,0,1,4])
            print u.T
            print U.getBlock()
        

    if ind == 1:
        print "Bond ID: Right Bond"
        U = get_unitary4bond_up(Ulist,L,L-2)
        u = np.einsum('ijkl,mknp->mijnpl',u13,u22).reshape(8,8)
        #print u13.reshape(4,4)
        #print u22.reshape(4,4)
        if dag is False:
            print "Numpy result: "
            print u
            print "My result: "
            print U.getBlock()
            print U
        else:
            U.permute([5,6,7,4,2,3],3)
            U.setLabel([8,9,10,4,0,1])
            print u.T
            print U.getBlock()

    if ind == 2:
        print "Bond ID: Odd Bond"
        U = get_unitary4bond_up(Ulist,L,1)
        u = np.einsum('ijkl,lmnp,qrms -> ijqrknps',u11,u21,u12).reshape(16,16)
        #print u11.reshape(4,4)
        #print u12.reshape(4,4)
        #print u21.reshape(4,4)

        if dag is False:
            print u
            print U.getBlock()
        else:
            U.permute([5,6,7,8,4,2,3,-4],4)
            U.setLabel([9,10,11,12,4,0,1,-4])
            print u.T
            print U.getBlock()
    


    if ind == 3:
        print "Bond ID: Even Bond"
        U = get_unitary4bond_up(Ulist,L,2)
        u = np.einsum('ijkl,nkmp,lqrs -> nijqmprs',u12,u21,u22).reshape(16,16)
        #print u12.reshape(4,4)
        #print u21.reshape(4,4)
        #print u22.reshape(4,4)
        if dag is False:
            print u
            print U.getBlock()
        else:
            U.permute([5,6,7,8,4,2,3,-4],4)
            U.setLabel([9,10,11,12,4,0,1,-4])
            print u.T
            print U.getBlock()
    



    # Test U-dagger construction
    
    

def test_UHU(ind):
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
    
    #Hm = get_Hbond_tensor(edge=False)
    #Hl = get_Hbond_tensor(edge=True, edgeside='L')
    #Hr = get_Hbond_tensor(edge=True, edgeside='R')
    
    
    hm = -1.* np.kron(sx,sx)-0.5*(np.kron(sz,si)+np.kron(si,sz))
    hl = -1.* np.kron(sx,sx)-0.5*np.kron(si,sz)-np.kron(sz,si)
    hr = -1.* np.kron(sx,sx)-np.kron(si,sz)-0.5*np.kron(sz,si)
    

    UHU = contract_UHU(6,Ulist)
    #print UHU[0] 
    
    
    if ind == 0: # LEFT EDGE
        u11_ = u11.reshape(4,4)
        uhu1 = np.dot(u11_.T,np.dot(hl,u11_))
        uhu1 = uhu1.reshape(2,2,2,2)
        u21T = (u21.reshape(4,4).T).reshape(2,2,2,2)
        uhu = np.einsum('ijkl,mknp,plrs->mijnrs',u21T,uhu1,u21).reshape(8,8)
        print uhu
        print UHU[0].getBlock()
    elif ind == 1: # Interior Odd Bonds
        hm_ = hm.reshape(2,2,2,2)
        u11T = (u11.reshape(4,4).T).reshape(2,2,2,2)
        u12T = (u12.reshape(4,4).T).reshape(2,2,2,2)
        u21T = (u21.reshape(4,4).T).reshape(2,2,2,2)
        uhu1 = np.einsum('ijkl,mnpq,lprs,kruv,sqwx->ijmnuvwx',u11T,u12T,hm_,u11,u12)
        uhu  = np.einsum('ijkl,mklnpqrs,qruv->mijnpuvs',u21T,uhu1,u21).reshape(16,16)
        print uhu
        print UHU[1].getBlock()

    elif ind == 2:
        u12_ = u12.reshape(4,4)
        uhu1 = np.dot(u12_.T,np.dot(hm, u12_)).reshape(2,2,2,2)
        u21T = (u21.reshape(4,4).T).reshape(2,2,2,2)
        u22T = (u22.reshape(4,4).T).reshape(2,2,2,2)
        uhu = np.einsum('ijkl,mnpq,lpxy,kxrs,yquv->ijmnrsuv',u21T,u22T,uhu1,u21,u22).reshape(16,16)
        print uhu
        print UHU[2].getBlock()

    else:
        u13_ = u13.reshape(4,4)
        uhu1 = np.dot(u13_.T,np.dot(hr,u13_)).reshape(2,2,2,2)
        u22T = (u22.reshape(4,4).T).reshape(2,2,2,2)
        uhu  = np.einsum('ijkl,lmqn,kqrs->ijmrsn',u22T,uhu1,u22).reshape(8,8)
        print uhu
        print UHU[-1].getBlock()
        
        
    
    

test_UHU(int(sys.argv[1]))
