import numpy as np
import pyUni10 as uni10



def print_banner(Model,L,L_lay,J,Fieldz,Method):

    print "#####################################################"
    print "#        Variational Unitary Tensor Network         #"
    print "#####################################################"
    print "#  Parameters                                       #"
    print "#  Model: %2d site "%L+Model+" "*(34-len(Model))+"#"
    print "#  Number of Layers: %1d"%(len(L_lay))+" "*30+"#"
    print "#  Couplings: J = %0.1f, Fieldz = %0.1f"%(J,Fieldz) + " "*17+"#"
    print "#  Optimization algorithm: %s"%Method+" "*(25-len(Method))+"#"
    print "#####################################################" 

######PAULI STRINGS######
def matSx():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [0.0, 1.0, 1.00, 0.0])
  return Mat 

def matSz():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [1.0, 0, 0, -1.0]);
  return Mat 

def matSy():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [0.0, -1.00, 1.00, 0.00]);
  return Mat 

def matIden():
    spin_t=0.5
    dimT = int(2*spin_t + 1)
    Mat=uni10.Matrix(dimT, dimT,[1,0,0,1])
    return Mat
        
def pauli2body():
    '''
    Tensors for 1-body Pauli strings. 
    '''
    bond_dim = 2
    iden = matIden()
    sx   = matSx()
    sy   = matSy()
    sz   = matSz()
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)
    pauli1b = [iden,sx,sy,sz]
    pauli2b = []
    for i in xrange(4):
      for j in xrange(4):
        mat = uni10.otimes(pauli1b[i],pauli1b[j])
        P = uni10.UniTensor([bdi,bdi,bdo,bdo])
        P.putBlock(mat)
        P.setLabel([0,1,2,3])
        pauli2b.append(P)
    return pauli2b


def pauli3body():
    '''
    Tensors for 3-body Pauli strings. 
    '''
    bond_dim = 2
    iden = matIden()
    sx   = matSx()
    sy   = matSy()
    sz   = matSz()
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)

    pauli1b = [iden,sx,sy,sz]
    pauli3b = []
    for i in xrange(4):
      for j in xrange(4):
        for k in xrange(4):
          mat = uni10.otimes(uni10.otimes(pauli1b[i],pauli1b[j]),pauli1b[k])
          P   = uni10.UniTensor([bdi,bdi,bdi,bdo,bdo,bdo])
          P.putBlock(mat)
          pauli3b.append(P)

    return pauli3b

def pauli4body():
    '''
    Tensors for 4-body Pauli strings. 
    '''
    bond_dim = 2
    iden = matIden()
    sx   = matSx()
    sy   = matSy()
    sz   = matSz()
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)

    pauli1b = [iden,sx,sy,sz]
    pauli4b = []
    pauli2b = []
    l2b = 16
    for i in xrange(4):
        for j in xrange(4):
            mat = uni10.otimes(pauli1b[i], pauli1b[j])
            pauli2b.append(mat)
    for i in xrange(l2b):
        for j in xrange(l2b):
            mat = uni10.otimes(pauli2b[i], pauli2b[j])
            P   = uni10.UniTensor([bdi,bdi,bdi,bdi,bdo,bdo,bdo,bdo])
            P.putBlock(mat)
            pauli4b.append(P)            

    return pauli4b

def pauli5body():
    '''
    Tensors for 5-body Pauli strings. 
    '''
    bond_dim = 2
    iden = matIden()
    sx   = matSx()
    sy   = matSy()
    sz   = matSz()
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)

    pauli1b = [iden,sx,sy,sz]
    pauli5b = []
    bonds = [bdi]*5 + [bdo]*5
    pauli2b = []
    l2b = 16
    for i in xrange(4):
      for j in xrange(4):
        mat = uni10.otimes(pauli1b[i],pauli1b[j])
        pauli2b.append(mat)
    for i in xrange(l2b):
      for j in xrange(l2b):
        for k in xrange(4):
          mat = uni10.otimes(uni10.otimes(pauli2b[i],pauli2b[j]), pauli1b[k])
          P   = uni10.UniTensor(bonds)
          P.putBlock(mat)
          pauli5b.append(P)

    return pauli5b

def pauli6body():
    '''
    Tensors for 6-body Pauli strings. 
    '''
    bond_dim = 2
    iden = matIden()
    sx   = matSx()
    sy   = matSy()
    sz   = matSz()
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)

    pauli1b = [iden,sx,sy,sz]
    pauli6b = []
    bonds = [bdi]*6 + [bdo]*6
    pauli3b = []
    
    for i in xrange(4):
      for j in xrange(4):
        for k in xrange(4):
            mat = uni10.otimes(uni10.otimes(pauli1b[i],pauli1b[j]),pauli1b[k])
            pauli3b.append(mat)
    l3b = 64
    for i in xrange(l3b):
        for j in xrange(l3b):
            mat = uni10.otimes(pauli3b[i],pauli3b[j])
            P   = uni10.UniTensor(bonds)
            P.putBlock(mat)
            pauli6b.append(P)

    return pauli6b

###### UNI10 to NUMPY ######

def uni2mat(uni,lmat):
    # uni10 matrix to numpy array
    try:
        M = uni.getBlock()
    except:
        M = uni
    npmat = np.zeros((lmat,lmat))
    for i in xrange(lmat):
        for j in xrange(lmat):
            npmat[i,j] = M[i*lmat+j]
    return npmat

def mat2uni(mat,lmat):
    mat = mat.reshape(lmat**2)
    mat = list(mat)
    unimat = uni10.Matrix(lmat,lmat,mat)
    return unimat

def vec2uni(v,L):
    bond_dim = 2
    vl = list(v)
    vl10 = uni10.Matrix(2**L,1,vl)
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)
    vtensor = uni10.UniTensor([bdi]*L,'GS')
    vtensor.putBlock(vl10)
    return vtensor

def get_uni10diag(uni,lmat):
    try:
        M = uni.getBlock()
    except:
        M = uni
    npdiag = np.zeros(lmat)
    for i in xrange(lmat):
        npdiag[i] = M[i*(lmat+1)]
    return npdiag


def str2int(mystr, base):
    l = len(mystr)
    str_n = mystr[::-1]
    myint = 0
    for i in xrange(l):
        myint += str_n[i]*(base**i)
    return int(myint)

###### ARITHMETIC ######
def diff_unitary(Ulist1, Ulist2):
    lsite = len(Ulist1)
    llayer = len(Ulist1[0])
    if (len(Ulist2) != lsite or len(Ulist2[0]) != llayer):
        raise Exception("The dimensions of the two unitaries do not match!")
    diff = 0.
    for i in xrange(lsite):
        for j in xrange(llayer):
            u1 = Ulist1[i][j].getBlock()
            u2 = Ulist2[i][j].getBlock()
            diff += (-1.*u1 + u2).norm()

    diff /= lsite*llayer
    return diff
            
    
    
