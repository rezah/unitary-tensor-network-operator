import pyUni10 as uni10
#import syss
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import random
import copy
import time
import mpo
import env
import optimize
################################## Set parameter #######################################
Model='Ising'               #could be: 'Ising' and or 'Heisenberg'
L=12                             #Numbers of particles should be even!
L_lay=[0,1]                  #Numbers of layers <= 5
d=2                              #pysical bond-dimension  
chi=3                            #bond-dimension of MPO

J=1.0                            #coupling, 1.0
Fieldz=1.0                            #Field in z direction, 0.6


W=8                              #random interval, [-W,W]
J2=0.0                           #coupling, 0.3
hz_list=[]                        # list of randomness 
U_delta=0.0                     #if it's zero, U_list is intialize by Identity
Method='CGploy'          #methods: CGarmjo,CGploy, SVD, SteepestDescent, SteepestDescentploy 
Randomness='Fixed'
Max_number_iteratoin_SVD=60        # maximum number of sweeps for SVD method 
Max_number_iteratoin_Steepest=60   # maximum number of sweeps for SteepestDescent method
Max_number_iteratoin_CG=60   # maximum number of sweeps for SteepestDescent method
Max_SVD_iteratoin=30               #maximum number of SVD iteration
Max_Steepest_iteratoin=30          #maximum number of SteepestDescent iteration
Max_CG_iteratoin=30          #maximum number of SteepestDescent iteration
######################################################
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


def Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list, Gamma):
 if Method is 'SVD':
  Length=Max_number_iteratoin_SVD
 elif Method is 'SteepestDescent':
  Length=Max_number_iteratoin_Steepest
 elif (Method is 'CGarmjo') or (Method is 'CGpoly') or (Method is 'SteepestDescentploy'):
  Length=Max_number_iteratoin_CG
 for j in xrange(Length):
  for i in xrange( L/2 ):
   #making Env right
   if i is 0:
    for m in xrange((L/2)-1):
     L_position1=(L/2)-1-m
     #print L_position
     Environment_Right[L_position1]=env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position1, d, Environment_Right)

   #print i
   L_position=i
   optimize.optimize_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list,Gamma)
   if i is not L/2 -1:
    Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)

#################################################################
def Initialize_function(L, L_lay):
 Env_Uni_inner=[]
 Gamma=[]
 for i in xrange(L/2):
  Env_Uni_inner.append([])
  Gamma.append([])
 for i in xrange(L/2):
  for j in xrange(len(L_lay)):
   Env_Uni_inner[i].append([None])
   Gamma[i].append(1.00)
 return Env_Uni_inner, Gamma

#####################################################################################
mpo.intialize_coupling(L, hz_list, W, Randomness, Model)

print "couplings", J, Fieldz, "L", L, "N_lay", len(L_lay)
######### retrun MPO list, 1: one-site MPO, 2: two-site MPO; mpo_uni10_bl, mpo_uni10_br, mpo_uni10_bl2, mpo_uni10_br2###
mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model )
trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L )
print trH2

U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta)

#####################UHU####################UHU####################UHU
## construct XX+hZ
#def get_Hbond_tensor(bond_dim=2, edge=False,edgeside='L'):
#
#    '''
#    XX + h Z
#    Result in a tensor:
#            _____
#       0 __|     |__ 2
#       1 __|  H  |__ 3
#           |_____| 
#
#    '''
#    # make matrix
#    sx   = matSx()
#    sy   = matSy()
#    sz   = matSz()
#    iden = matIden()
#    Hbond_mat  = (-1.*J)*uni10.otimes(sx,sx)
#    if edge:
#        if edgeside is 'L':
#            Hbond_mat += (-1.0*Fieldz)*uni10.otimes(sz,iden) 
#            Hbond_mat += (-0.5*Fieldz)*uni10.otimes(iden,sz)
#        elif edgeside is 'R':
#            Hbond_mat += (-0.5*Fieldz)*uni10.otimes(sz,iden) 
#            Hbond_mat += (-1.0*Fieldz)*uni10.otimes(iden,sz)
#    else:
#        Hbond_mat += (-0.5*Fieldz)*uni10.otimes(sz,iden) 
#        Hbond_mat += (-0.5*Fieldz)*uni10.otimes(iden,sz)
#    # construct tensor
#    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
#    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)
#
#    Hbond = uni10.UniTensor([bdi,bdi,bdo,bdo],'Hbond')
#    Hbond.putBlock(Hbond_mat)
#    
#    return Hbond
#
#'''
##test1
#H1 = get_Hbond_tensor()
#H2 = get_Hbond_tensor(edge=True)
#
#print "H1:"
#print H1
#H2.setLabel([2,3,4,5])
#print "H2:"
#print H2
#print "+++++++"
#print H1*H2
#H3 = copy.copy(H2)
#H3.permute([4,5,2,3],2)
#H3.setLabel([2,3,4,5])
#print "H3"
#print H3
#print "=========="
#print H1*H3
#exit()
#'''         
#
#
#def get_unitary4bond_up(Ulist,bond_idx,bond_dim=2,edge=False,edgeside='L'):
#
#    # TODO Use bond_idx to detect edge bonds, how to detect right edge?
#    # Ulist is a L/2 x 2 list, 2 is the number of layers 
#    l_layer = len(Ulist[0])
#    if l_layer != 2:
#        raise ValueError("Please make sure your unitary list has TWO layers!")
#    
#    if edge:
#        if edgeside is 'L':
#            '''
#               _____
#          2 __|     |__ 5
#          3 __|  U  |__ 6
#          4 __|     |__ 7 
#              |_____| 
#        
#           Where 2,3 will combine with the Hamiltonian,
#           4 will combine with U_down
#            '''
#            U0 = copy.copy(Ulist[0][0])
#            U1 = copy.copy(Ulist[0][1])
#            U0.setLabel([0,1,2,3])
#            U1.setLabel([3,4,5,6])
#            U = U0*U1
#            U.setLabel([2,3,4,5,6,7])
#
#        elif edgeside is 'R':
#            '''
#               _____
#          4 __|     |__ 5
#          2 __|  U  |__ 6
#          3 __|     |__ 7 
#              |_____| 
#        
#           Where 2,3 will combine with the Hamiltonian,
#           4 will combine with U_down
#            '''
#
#            U0 = copy.copy(Ulist[-1][0])
#            U1 = copy.copy(Ulist[-1][1])
#            U0.setLabel([0,1,2,3])
#            U1.setLabel([4,2,5,6])
#            U = U0*U1
#            U.setLabel([4,2,3,5,6,7])
#        else:
#            raise ValueError("The edgeside argument can only be 'L' or 'R'! ")
#    else:
#
#        '''
#            _____
#       4 __|     |__ 5
#       2 __|     |__ 6
#       3 __|  U  |__ 7 
#      -4 __|     |__ 8        
#           |_____| 
#     
#        Where 2,3 will combine with the Hamiltonian,
#        4,-4 will combine with U_down
#        '''
#
#        if(bond_idx//2==0):
#            l = bond_idx/2
#            U0  = copy.copy(Ulist[l][0])
#            U1L = copy.copy(Ulist[l-1][1])
#            U1R = copy.copy(Ulist[l][1])
#            U0.setLabel([0,1,2,3])
#            U1L.setLabel([3,4,5,6])
#            U1R.setLabel([7,2,8,9])
#            U = (U0*U1L)*U1R
#            U.setLabel([4,2,3,-4,5,6,7,8]) 
#
#        else:
#            l = bond_idx//2
#            U0L = copy.copy(Ulist[l][0])
#            U0R = copy.copy(Ulist[l+1][0])
#            U1  = copy.copy(Ulist[l][1])
#            U0L.setLabel([0,1,2,3])
#            U0R.setLabel([4,5,6,7])
#            U1.setLabel([7,2,8,9])
#            U = U0L*(U0R*U1)
#            U.setLabel([4,2,3,-4,5,6,7,8]) 
#    return U
#
#'''
##test2
#U = get_unitary4bond_up(U_list,0,edge=True) 
#print U
#'''
#
#def contract_UHU(L, Ulist):
#    '''
#    Return a list of U^+HU.
#    L is the lattice length.
#    '''
#
#    UHU = []
#    # left edge
#    Hbond = get_Hbond_tensor(edge=True,edgeside='L') 
#    U_up  = get_unitary4bond_up(Ulist,0,edge=True,edgeside='L')
#    U_dn  = copy.copy(U_up)
#    U_dn.permute([5,6,7,2,3,4],3)
#    U_dn.setLabel([8,9,10,0,1,4])
#    Hbond_n = (U_dn*Hbond)*U_up
#    UHU.append(Hbond_n)
#
#    for bidx in range(1,L-2):
#        Hbond = get_Hbond_tensor(edge=False)
#        U_up  = get_unitary4bond_up(Ulist,bidx)
#        U_dn  = copy.copy(U_up)
#        U_dn.permute([5,6,7,8,4,2,3,-4],4)
#        U_dn.setLabel([9,10,11,12,4,0,1,-4])
#        Hbond_n = (U_dn*Hbond)*U_up
#        UHU.append(Hbond_n)
#
#    Hbond = get_Hbond_tensor(edge=True,edgeside='R')
#    U_up  = get_unitary4bond_up(Ulist,L-1,edge=True,edgeside='R')
#    U_dn  = copy.copy(U_up)
#    U_dn.permute([5,6,7,4,2,3],3)
#    U_dn.setLabel([8,9,10,4,0,1])
#    Hbond_n = (U_dn*Hbond)*U_up
#    UHU.append(Hbond_n)
#
#    return UHU
#    

import qsampl
UHU = qsampl.contract_UHU(4,U_list)
print UHU
    
            
        

            
            
     
            

    



            







    
