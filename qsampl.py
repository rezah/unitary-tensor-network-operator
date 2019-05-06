import pyUni10 as uni10
import copy
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
    XX + h Z
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

    for bidx in range(1,L-2):
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
def pauli2body():
    # TODO figure out an easier way to print the labels
    '''
    Tensors for 2-body Pauli strings. 
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
    pauli2b_labels = [['II','IX','IY','IZ'],\
                      ['XI','XX','XY','XZ'],\
                      ['YI','YX','YY','YZ'],\
                      ['ZI','ZX','ZY','ZZ']]            
    for i in range (4):
        for j in range(4):
            mat = uni10.otimes(pauli1b[i],pauli1b[j])
            P   = uni10.UniTensor([bdi,bdi,bdo,bdo],pauli2b_labels[i][j])
            P.putBlock(mat)
            pauli2b.append(P) 

    return pauli2b, pauli2b_labels

def pauli_table():
    '''
    Construct the Pauli multiplication table.
    Return a 4x4 signed Matrix. 
    +/- index.
    '''
    table = np.zeros((4,4),dtype=np.int32)
    for i in range(4):
        table[i,i] = 0
        table[0,i] = i
    table[1,2] = 3
    table[1,3] = 2
    table[2,3] = 1
    table = table + table.T

    table[2,1] *= -1
    table[3,1] *= -1
    table[3,2] *= -1
    
    return table
            
    
def unitary2pauli(U,P,lp):
    '''
    U  - unitary
    P  - Pauli string
    lp - length of Pauli string 
    '''
    U_ = copy.copy(U)
    U_.setLabel([0,1,2,3])
    trUP = (U_*P).getBlock().sum()/(2.*lp)
    return trUP

def ulist2pauli(Ulist,lp):
    '''
    Turn the Ulist into a list of vectors.
    '''
    #TODO Call pauli2body() in a higher fuction to reduce the number of calls
    ptensor, plabel = pauli2body() 
    l0 = len(Ulist)
    l1 = len(Ulist[0])
    np2 = len(ptensor)
    print np2
    
    Ucoef = []
    
    if l1 != 2:
        raise ValueError("Only 2-layer U_list is accepted!") 
    for i in range(l0):
        ucoef_ = []
        for j in range(l1):
            coef = np.zeros(np2)
            U = Ulist[i][j]
            for k in range(np2):
                coef[k] = unitary2pauli(U,ptensor[k],lp)
                if (abs(coef[k])<1e-10):
                    coef[k] = 0.   # setting small coefficients to zero
            ucoef_.append(coef)
        Ucoef.append(ucoef_)     

    return Ucoef



            
            
            
   


    
    
    
    

            
        

            
            
 
