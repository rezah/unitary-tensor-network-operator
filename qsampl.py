import pyUni10 as uni10
import copy

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
    sy   = matSy()
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

'''
#test1
H1 = get_Hbond_tensor()
H2 = get_Hbond_tensor(edge=True)

print "H1:"
print H1
H2.setLabel([2,3,4,5])
print "H2:"
print H2
print "+++++++"
print H1*H2
H3 = copy.copy(H2)
H3.permute([4,5,2,3],2)
H3.setLabel([2,3,4,5])
print "H3"
print H3
print "=========="
print H1*H3
exit()
'''         


def get_unitary4bond_up(Ulist,bond_idx,bond_dim=2,edge=False,edgeside='L'):

    # TODO Use bond_idx to detect edge bonds, how to detect right edge?
    # Ulist is a L/2 x 2 list, 2 is the number of layers 
    l_layer = len(Ulist[0])
    if l_layer != 2:
        raise ValueError("Please make sure your unitary list has TWO layers!")
    
    if edge:
        if edgeside is 'L':
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
            U.setLabel([2,3,4,5,6,7])

        elif edgeside is 'R':
            '''
               _____
          4 __|     |__ 5
          2 __|  U  |__ 6
          3 __|     |__ 7 
              |_____| 
        
           Where 2,3 will combine with the Hamiltonian,
           4 will combine with U_down
            '''

            U0 = copy.copy(Ulist[-1][0])
            U1 = copy.copy(Ulist[-1][1])
            U0.setLabel([0,1,2,3])
            U1.setLabel([4,2,5,6])
            U = U0*U1
            U.setLabel([4,2,3,5,6,7])
        else:
            raise ValueError("The edgeside argument can only be 'L' or 'R'! ")
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

        if(bond_idx//2==0):
            l = bond_idx/2
            U0  = copy.copy(Ulist[l][0])
            U1L = copy.copy(Ulist[l-1][1])
            U1R = copy.copy(Ulist[l][1])
            U0.setLabel([0,1,2,3])
            U1L.setLabel([3,4,5,6])
            U1R.setLabel([7,2,8,9])
            U = (U0*U1L)*U1R
            U.setLabel([4,2,3,-4,5,6,7,8]) 

        else:
            l = bond_idx//2
            U0L = copy.copy(Ulist[l][0])
            U0R = copy.copy(Ulist[l+1][0])
            U1  = copy.copy(Ulist[l][1])
            U0L.setLabel([0,1,2,3])
            U0R.setLabel([4,5,6,7])
            U1.setLabel([7,2,8,9])
            U = U0L*(U0R*U1)
            U.setLabel([4,2,3,-4,5,6,7,8]) 
    return U

'''
#test2
U = get_unitary4bond_up(U_list,0,edge=True) 
print U
'''

def contract_UHU(L, Ulist):
    '''
    Return a list of U^+HU.
    L is the lattice length.
    '''

    UHU = []
    # left edge
    Hbond = get_Hbond_tensor(edge=True,edgeside='L') 
    U_up  = get_unitary4bond_up(Ulist,0,edge=True,edgeside='L')
    U_dn  = copy.copy(U_up)
    U_dn.permute([5,6,7,2,3,4],3)
    U_dn.setLabel([8,9,10,0,1,4])
    Hbond_n = (U_dn*Hbond)*U_up
    UHU.append(Hbond_n)

    for bidx in range(1,L-2):
        Hbond = get_Hbond_tensor(edge=False)
        U_up  = get_unitary4bond_up(Ulist,bidx)
        U_dn  = copy.copy(U_up)
        U_dn.permute([5,6,7,8,4,2,3,-4],4)
        U_dn.setLabel([9,10,11,12,4,0,1,-4])
        Hbond_n = (U_dn*Hbond)*U_up
        UHU.append(Hbond_n)

    Hbond = get_Hbond_tensor(edge=True,edgeside='R')
    U_up  = get_unitary4bond_up(Ulist,L-1,edge=True,edgeside='R')
    U_dn  = copy.copy(U_up)
    U_dn.permute([5,6,7,4,2,3],3)
    U_dn.setLabel([8,9,10,4,0,1])
    Hbond_n = (U_dn*Hbond)*U_up
    UHU.append(Hbond_n)

    return UHU
    
    
            
        

            
            
 
