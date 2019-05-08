import pyUni10 as uni10
import copy
import itertools
import numpy as np

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

####################UHU####################UHU####################UHU
def get_Hbond_tensor(J=1.0,Fieldz=1.0, bond_dim=2, edge=False,edgeside='L'):

    '''
    -J XX - h Z    (h = Fieldz)
    Result in a tensor:
            _____
       0 __|     |__ 2
       1 __|  H  |__ 3
           |_____| 

    '''
    # make matrix
    sx   = matSx()
    #sy   = matSy()
    sz   = matSz()
    iden = matIden()
    Hbond_mat  = (-1.*J)*uni10.otimes(sx,sx)
    if edge:
        if edgeside is 'L':
            Hbond_mat += (-1.0*Fieldz)*uni10.otimes(sz,iden) 
            Hbond_mat += (-0.5*Fieldz)*uni10.otimes(iden,sz)
        elif edgeside is 'R':
            Hbond_mat += (-0.5*Fieldz)*uni10.otimes(sz,iden) 
            Hbond_mat += (-1.0*Fieldz)*uni10.otimes(iden,sz)
    else:
        Hbond_mat += (-0.5*Fieldz)*uni10.otimes(sz,iden) 
        Hbond_mat += (-0.5*Fieldz)*uni10.otimes(iden,sz)
    # construct tensor
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)

    Hbond = uni10.UniTensor([bdi,bdi,bdo,bdo],'Hbond')
    Hbond.putBlock(Hbond_mat)
    
    return Hbond


def get_unitary4bond_up(Ulist,L,bond_idx,bond_dim=2):


    #TODO Combine some tensors in advance and save them to speed up
    # Ulist is a L/2 x 2 list, 2 is the number of layers 
    # bond_idx ranges from 0 to L-2
    NU1 = L/2
    l_layer = len(Ulist[0])
    if l_layer != 2:
        raise ValueError("Please make sure your unitary list has TWO layers!")
    
    if (bond_idx == 0):   # Left Edge
        '''
              _____
         2 __|     |__ 5
         3 __|  U  |__ 6
         4 __|     |__ 7 
             |_____| 
    
        Where 2,3 will combine with the Hamiltonian,
        4 will combine with U_down
        '''
        U0 = copy.copy(Ulist[0][0])
        U1 = copy.copy(Ulist[0][1])
        U0.setLabel([0,1,2,3])
        U1.setLabel([3,4,5,6])
        U = U0*U1
        U.permute([0,1,4,2,5,6],3)
        U.setLabel([2,3,4,5,6,7])

    elif (bond_idx == L-2):   # Right Edge
        '''
              _____
         4 __|     |__ 5
         2 __|  U  |__ 6
         3 __|     |__ 7 
             |_____| 
       
        Where 2,3 will combine with the Hamiltonian,
        4 will combine with U_down
        '''

        U0 = copy.copy(Ulist[NU1-1][0])
        U1 = copy.copy(Ulist[NU1-2][1])
        U0.setLabel([0,1,2,3])
        U1.setLabel([4,2,5,6])
        U = U0*U1
        U.permute([4,0,1,5,6,3],3)
        U.setLabel([4,2,3,5,6,7])
    else:


        '''
            _____
       4 __|     |__ 5
       2 __|     |__ 6
       3 __|  U  |__ 7 
      -4 __|     |__ 8        
           |_____| 
     
        Where 2,3 will combine with the Hamiltonian,
        4,-4 will combine with U_down
        '''

        if(bond_idx%2 == 0):   # Even Bond
            l = bond_idx/2
            U0  = copy.copy(Ulist[l][0])
            U1L = copy.copy(Ulist[l-1][1])
            U1R = copy.copy(Ulist[l][1])
            U0.setLabel([0,1,2,3])
            U1L.setLabel([4,2,5,6])
            U1R.setLabel([3,7,8,9])
            U = (U0*U1L)*U1R
            U.permute([4,0,1,7,5,6,8,9],4)
            U.setLabel([4,2,3,-4,5,6,7,8]) 

        else:                # Odd Bond
            l = bond_idx//2
            U0L = copy.copy(Ulist[l][0])
            U0R = copy.copy(Ulist[l+1][0])
            U1  = copy.copy(Ulist[l][1])
            U0L.setLabel([0,1,2,3])
            U0R.setLabel([4,5,6,7])
            U1.setLabel([3,6,8,9])
            U = U0L*(U0R*U1)
            U.permute([0,1,4,5,2,8,9,7],4)
            U.setLabel([4,2,3,-4,5,6,7,8]) 

    return U

def contract_UHU(L, Ulist):
    '''
    Return a list of U^+HU.
    L is the lattice length.
    '''
    #TODO Consider when U has imaginary parts

    UHU = []
    # left edge
    Hbond = get_Hbond_tensor(edge=True,edgeside='L') 
    U_up  = get_unitary4bond_up(Ulist,L,0)
    U_dn  = copy.copy(U_up)
    U_dn.permute([5,6,7,2,3,4],3)
    U_dn.setLabel([8,9,10,0,1,4])
    Hbond_n = (U_dn*Hbond)*U_up
    Hbond_n.setLabel([0,1,2,3,4,5])
    UHU.append(Hbond_n)

    for bidx in xrange(1,L-2):
        Hbond = get_Hbond_tensor(edge=False)
        U_up  = get_unitary4bond_up(Ulist,L,bidx)
        U_dn  = copy.copy(U_up)
        U_dn.permute([5,6,7,8,4,2,3,-4],4)
        U_dn.setLabel([9,10,11,12,4,0,1,-4])
        Hbond_n = (U_dn*Hbond)*U_up
        Hbond_n.setLabel([0,1,2,3,4,5,6,7])
        UHU.append(Hbond_n)

    Hbond = get_Hbond_tensor(edge=True,edgeside='R')
    U_up  = get_unitary4bond_up(Ulist,L,L-2)
    U_dn  = copy.copy(U_up)
    U_dn.permute([5,6,7,4,2,3],3)
    U_dn.setLabel([8,9,10,4,0,1])
    Hbond_n = (U_dn*Hbond)*U_up
    Hbond_n.setLabel([0,1,2,3,4,5])
    UHU.append(Hbond_n)

    return UHU

       
######PAULI######PAULI######PAULI######PAULI######PAULI######PAULI
def trUP(U,P,lp):
    '''
    U  - unitary
    P  - Pauli string
    lp - length of Pauli string 
    '''
    #U_ = copy.copy(U)
    #U_.setLabel([0,1,2,3])
    tr = (U*P).getBlock().sum()/(2.**lp)
    return tr

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
    for i in xrange(4):
      for j in xrange(4):
        for k in xrange(4):
          for l in xrange(4):
            mat = uni10.otimes(uni10.otimes(uni10.otimes(pauli1b[i],pauli1b[j]),pauli1b[k]),pauli1b[l])
            P   = uni10.UniTensor([bdi,bdi,bdi,bdi,bdo,bdo,bdo,bdo])
            P.putBlock(mat)
            pauli4b.append(P)

    return pauli4b


def UHU2pauli(L,Ulist):
    # TODO make very small terms zero so that no measurement is needed
    #      can be done when measuring too
    UHU = contract_UHU(L, Ulist) 
    pauli3b = pauli3body()
    pauli4b = pauli4body()
    l3b = 64
    l4b = 256
    pcoefs = []

    coef3 = np.zeros(l3b)
    for i in xrange(l3b):
        coef3[i] = trUP(UHU[0],pauli3b[i],3)
    pcoefs.append(coef3)
    for l in xrange(1,L-2):
        coef4 = np.zeros(l4b)
        for i in xrange(l4b):
            coef4[i] = trUP(UHU[l],pauli4b[i],4)
        pcoefs.append(coef4)

    coef3 = np.zeros(l3b)
    for i in xrange(l3b):
        coef3[i] = trUP(UHU[L-2],pauli3b[i],3) 
    pcoefs.append(coef3)

    return pcoefs


##########GROUND#STATE##########GROUND#STATE##########GROUND#STATE
def vector2uni10(v,L):
    bond_dim = 2
    vl = list(v)
    vl10 = uni10.Matrix(2**L,1,vl)
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)
    vtensor = uni10.UniTensor([bdi]*L,'GS')
    vtensor.putBlock(vl10)
    return vtensor

 
def ising_ED(L,J=1.0,Fieldz=1.0):
    '''
    -J XX - h Z    (h = Fieldz)
    '''

    if (L > 12):
        raise ValueError("Exact diagonalization can only handle up to 12 sites!")
    if (L%2 != 0):
        raise ValueError("Only even number of sites is accepted!")

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
    gs = vector2uni10(ev[:,0],L)
    return ew[0], ev[:,0] #TODO change ev[:,0] to gs after test

##############Upsi##############Upsi##############Upsi##############Upsi
def UdagPsi(Ulist, psi, L):
    lu = L/2
    if len(Ulist[0]) != 2:
        raise ValueError("Only 2-layer unitray towers are accepted!")
    v = copy.copy(psi)
    # Applying the first layer
    label = []
    for i in range(lu):
        s = i*2
        U0 = copy.copy(Ulist[i][0])
        U0.permute([2,3,0,1],2)
        U0.setLabel([s+L,s+L+1,s,s+1])
        v = U0*v
        label.append([s,s+1])
    label = label[::-1]
    label = list(itertools.chain(*label))
    v.setLabel(label)
    v.permute(xrange(L),L)

    # Applying the second layer
    label = []
    for i in range(lu-1):
        s=i*2+1
        U1 = copy.copy(Ulist[i][1])
        U1.permute([2,3,0,1],2)
        U1.setLabel([s+L,s+L+1,s,s+1])
        v = U1*v
        label.append([s,s+1])
    label = label[::-1]
    label = list(itertools.chain(*label))
    label = label + [0,L-1]
    v.setLabel(label)
    v.permute(xrange(L),L)
    
    return v
        
    


        
    

    
     




   




    
     
    

    
    




