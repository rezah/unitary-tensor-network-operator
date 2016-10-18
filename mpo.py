import pyUni10 as uni10
#import sys
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import random 
import copy
import env


##################### Spin 1/2 Pauli Operators ###############################################
def matS0():
  return uni10.Matrix(2, 2, [1.0, 0.0, 0.0, 1.0])

def matSx():
  return uni10.Matrix(2, 2, [0.0, 0.50, 0.50, 0.0])

def matSy():
  return uni10.Matrix(2, 2, [0.0, -0.50, 0.50, 0.0])

def matSz():
  return uni10.Matrix(2, 2, [0.50, 0.0, 0.0, -0.50])
#############################################################################3



def intialize_coupling(L, hz_list, W,Randomness): 
 if Randomness is 'Fixed':
  random.seed(3)
 for i in xrange(L):
  hz_list.append(random.uniform(-W,W))


def print_trH2(L, hz_list,J): 
 hz_list_2=[x**2 for x in hz_list]
 trH2=(2.00**L)*( J**2 * 3.0/16.00 * (L-1) + 1.0/4.0 * sum(hz_list_2))
 variance_Energy=(2.00**L)*( J**2 /16.00 * (L-1) + 1.0/4.0 * sum(hz_list_2))
 return trH2
 print trH2


def avarage_Energy_power_2(L, hz_list,J): 
 hz_list_2=[x**2 for x in hz_list]
 trH2=(2.00**L)*( J**2 * 3.0/16.00 * (L-1) + 1.0/4.0 * sum(hz_list_2))
 variance_Energy=(2.00**L)*( J**2 /16.00 * (L-1) + 1.0/4.0 * sum(hz_list_2))
 return variance_Energy
 

def variance_Energy_function(L, hz_list,J,Avarage_E_power2): 
 hz_list_2=[x**2 for x in hz_list]
 trH2=(2.00**L)*( J**2 * 3.0/16.00 * (L-1) + 1.0/4.0 * sum(hz_list_2))
 print 'variance_Energy=%d - %d' %(trH2, Avarage_E_power2)  
 return (trH2 - Avarage_E_power2)/(2.00**L)


def Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left):
 mpo_U_list_up= make_mpo_U_list(U_list, L_lay, L, 'up')
 mpo_U_list_up= make_mpo_U_list(U_list, L_lay, L, 'up')
 mpo_U_list_down= make_mpo_U_list(U_list, L_lay, L, 'down')
 Environment_Left=[None]*(L/2)
 Environment_Right=[None]*(L/2)

 for i in xrange((L/2)-1):
  L_position=i
  Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)

 for i in xrange((L/2)-1):
  L_position=(L/2)-1-i
  #print L_position
  Environment_Right[L_position]=env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Right)

 return (Environment_Right[1] * Environment_Left[0])[0] 






def make_mpo_H(L, J, hz_list ):
 d=2
 chi=5
 s0=matS0()
 sx=matSx()
 sy=matSy()
 sz=matSz()


###############  Define Boundry of MPO   ################################3
 mpo_matrix_bl=uni10.Matrix(1,chi)
 mpo_matrix_br=uni10.Matrix(chi,1)
 bdi_mpo = uni10.Bond(uni10.BD_IN, chi);
 bdo_mpo = uni10.Bond(uni10.BD_OUT, chi);
 mpo_uni10_bl=uni10.UniTensor([bdo_mpo], "mpo_uni10_bl")
 mpo_uni10_br=uni10.UniTensor([bdi_mpo], "mpo_uni10_br")


 for i in xrange(1):
   for j in xrange(chi):
    mpo_matrix_bl[i*chi+j]=0
    if(j==4 and i==0):
     mpo_matrix_bl[i*chi+j]=1


 for i in xrange(chi):
   for j in xrange(1):
    mpo_matrix_br[i*1+j]=0
    if(i==0 and j==0):
     mpo_matrix_br[i*1+j]=1

 mpo_uni10_bl.putBlock(mpo_matrix_bl)
 mpo_uni10_br.putBlock(mpo_matrix_br)
 
 mpo_uni10_bl_prime=copy.copy(mpo_uni10_bl)
 mpo_uni10_br_prime=copy.copy(mpo_uni10_br)

 mpo_uni10_bl_prime.setLabel([1])
 mpo_uni10_br_prime.setLabel([1])

 mpo_uni10_bl_2=uni10.contract(mpo_uni10_bl,mpo_uni10_bl_prime)
 mpo_uni10_br_2=uni10.contract(mpo_uni10_br,mpo_uni10_br_prime)

 mpo_uni10_bl_2.permute([0,1],0)
 mpo_uni10_br_2.permute([0,1],2)

 mpo_uni10_bl_2.combineBond([0,1])
 mpo_uni10_br_2.combineBond([0,1])
 mpo_uni10_bl_2.setName('mpo_uni10_bl_2')
 mpo_uni10_br_2.setName('mpo_uni10_br_2')
##############################################

 mpo_uni10_list_boundry=[]
 mpo_uni10_list_boundry.append(mpo_uni10_bl)
 mpo_uni10_list_boundry.append(mpo_uni10_br)
 mpo_uni10_list_boundry.append(mpo_uni10_bl_2)
 mpo_uni10_list_boundry.append(mpo_uni10_br_2)

#################### Define MPO ########################
 mpo_matrix=uni10.Matrix(chi*d,chi*d)


 bdi_spin = uni10.Bond(uni10.BD_IN, d);
 bdo_spin = uni10.Bond(uni10.BD_OUT, d);

 mpo_uni10_1_list=[]
 mpo_uni10_4_list=[]
 
 for num in xrange(L):
  mpo_matrix.set_zero()
  for mu in xrange(chi):
   for mup in xrange(chi):
    for s in xrange(d):
     for sp in xrange(d):
        if(mu==0 and mup==0  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=s0[s*d+sp]
        if(mu==1 and mup==0  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=sx[s*d+sp]
        if(mu==2 and mup==0  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=sy[s*d+sp]
        if(mu==3 and mup==0  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=sz[s*d+sp]
        if(mu==4 and mup==0  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=float(hz_list[num])*sz[s*d+sp]
        if(mu==4 and mup==1  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=float(J)*sx[s*d+sp]
        if(mu==4 and mup==2  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=-float(J)*sy[s*d+sp]
        if(mu==4 and mup==3  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=float(J)*sz[s*d+sp]
        if(mu==4 and mup==4  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=s0[s*d+sp]

 
  mpo_uni10=uni10.UniTensor([bdi_mpo, bdi_spin, bdo_mpo,bdo_spin], "mpo_uni10_1")
  
  mpo_uni10.putBlock(mpo_matrix)
  mpo_uni10_1_list.append(mpo_uni10)
  
 mpo_uni10_1_list_prime=[copy.copy(mpo_uni10_1_list[i]) for i in xrange(len(mpo_uni10_1_list)) ]

 for i in xrange(0,L,2):
   mpo_uni10_1_list_prime[i+1].setLabel([2,4,5,6])
   mpo_uni10_cotraction=uni10.contract(mpo_uni10_1_list_prime[i],mpo_uni10_1_list_prime[i+1] )
   mpo_uni10_cotraction.permute([0,1,4,5,3,6],3)
   mpo_uni10_cotraction.combineBond([1,4])
   mpo_uni10_cotraction.combineBond([3,6])
   mpo_uni10_cotraction.setLabel([0,1,2,3])
   mpo_uni10_cotraction.setName('mpo_uni10_4')
   mpo_uni10_4_list.append(mpo_uni10_cotraction)
   
   
 ##printmpo_uni10_1_list[1].printDiagram()
 
 return mpo_uni10_1_list, mpo_uni10_4_list, mpo_uni10_list_boundry



def contraction_MPO_trH2(mpo_uni10_1_list, mpo_uni10_4_list,  mpo_boundy_list):
 mpo_uni10_2_list=[]
 mpo_uni10_3_list=[]
 mpo_uni10_5_list=[]
 mpo_uni10_6_list=[]
 mpo_uni10_1_listN=[copy.copy(mpo_uni10_1_list[i])  for i in xrange(len(mpo_uni10_1_list)) ] 
 mpo_uni10_4_listN=[copy.copy(mpo_uni10_4_list[i]) for i in xrange(len(mpo_uni10_4_list)) ]
 L=len(mpo_uni10_1_listN)
 for num in xrange(L):
  mpo_uni10_prime=copy.copy(mpo_uni10_1_listN[num])
  mpo_uni10_prime.setLabel([4,5,6,1])
  mpo_uni10_2=uni10.contract(mpo_uni10_1_listN[num],mpo_uni10_prime)

  mpo_uni10_2.permute([0,4,5,2,6,3],3)
  mpo_uni10_2.combineBond([0,4])
  mpo_uni10_2.combineBond([2,6])
  mpo_uni10_2.setLabel([0,1,2,3])

  mpo_uni10_2.setName('mpo_uni10_2')
  mpo_uni10_2_list.append(mpo_uni10_2)
  
  mpo_uni10_2_prime=copy.copy(mpo_uni10_2)
  mpo_uni10_2_prime.partialTrace(1,3)
  mpo_uni10_2_prime.permute([0,2],1)
  mpo_uni10_2_prime.setLabel([0,1])
  mpo_uni10_2_prime.setName("mpo_uni10_3")
  mpo_uni10_3_list.append(mpo_uni10_2_prime)

 ##printmpo_uni10_list[0], mpo_uni10_list[1]
 ##printmpo_uni10_2_list[0], mpo_uni10_2_list[1]
 ##printmpo_uni10_3_list[0], mpo_uni10_3_list[1]

 mpo_uni10_3_list_prime=[copy.copy(mpo_uni10_3_list[i])  for i in xrange(len(mpo_uni10_3_list)) ]

 for i in xrange(L-1):
   if(i==0):
    mpo_uni10_3_list_prime[i].setLabel([0,1])
    mpo_uni10_3_list_prime[i+1].setLabel([1,2])
    mpo_uni10_ctraction=uni10.contract(mpo_uni10_3_list_prime[i], mpo_uni10_3_list_prime[i+1])
    mpo_uni10_ctraction.permute([0,2],1)
    mpo_uni10_ctraction.setLabel([0,1])
    mpo_uni10_ctraction.setName('mpo_uni10_ctraction')
    mpo_uni10_3_list_prime[0]=mpo_uni10_ctraction
   else:
    mpo_uni10_3_list_prime[i+1].setLabel([1,2])
    mpo_uni10_ctraction=uni10.contract(mpo_uni10_3_list_prime[0], mpo_uni10_3_list_prime[i+1])
    mpo_uni10_ctraction.permute([0,2],1)
    mpo_uni10_ctraction.setLabel([0,1])
    mpo_uni10_ctraction.setName('mpo_uni10_ctraction')
    mpo_uni10_3_list_prime[0]=mpo_uni10_ctraction

    

 #############  Testing calcualtion of trace H^{2} ###########################
 mpo_uni10_br_2=copy.copy(mpo_boundy_list[3])
 mpo_uni10_br_2.setLabel([1])
 trH2=mpo_boundy_list[2]*mpo_uni10_3_list_prime[0]*mpo_uni10_br_2
 #printtrH2[0]
 ###########################################################################



 mpo_uni10_1_list_prime=mpo_uni10_1_list[:]

 for i in xrange(0,L,2):
  
  mpo_uni10_cotraction_prime=copy.copy(mpo_uni10_4_listN[i/2])
  mpo_uni10_cotraction_prime.setLabel([4,5,6,1])

  mpo_uni10_contraction2=uni10.contract(mpo_uni10_4_listN[i/2], mpo_uni10_cotraction_prime)
  
  mpo_uni10_contraction2.permute([0,4,5,2,6,3],3)
  mpo_uni10_contraction2.combineBond([0,4])
  mpo_uni10_contraction2.combineBond([2,6])
  mpo_uni10_contraction2.setLabel([0,1,2,3])

  mpo_uni10_contraction2.setName('mpo_uni10_5')
  mpo_uni10_5_list.append(mpo_uni10_contraction2)
    
  mpo_uni10_cotraction3=copy.copy(mpo_uni10_contraction2)
  mpo_uni10_cotraction3.partialTrace(1,3)
  mpo_uni10_cotraction3.permute([0,2],1)
  mpo_uni10_cotraction3.setLabel([0,1])
  mpo_uni10_cotraction3.setName('mpo_uni10_6')
  mpo_uni10_6_list.append(mpo_uni10_cotraction3)
  ##printmpo_uni10_5_list[0], mpo_uni10_6_list[0] 
  



 mpo_uni10_6_list_prime=mpo_uni10_6_list[:]

 for i in xrange(len(mpo_uni10_6_list_prime)-1):
   if(i==0):
    mpo_uni10_6_list_prime[i].setLabel([0,1])
    mpo_uni10_6_list_prime[i+1].setLabel([1,2])
    mpo_uni10_ctraction=uni10.contract(mpo_uni10_6_list_prime[i], mpo_uni10_6_list_prime[i+1])
    mpo_uni10_ctraction.permute([0,2],1)
    mpo_uni10_ctraction.setLabel([0,1])
    mpo_uni10_ctraction.setName('mpo_uni10_ctraction')
    mpo_uni10_6_list_prime[0]=mpo_uni10_ctraction
   else:
    mpo_uni10_6_list_prime[i+1].setLabel([1,2])
    mpo_uni10_ctraction=uni10.contract(mpo_uni10_6_list_prime[0], mpo_uni10_6_list_prime[i+1])
    mpo_uni10_ctraction.permute([0,2],1)
    mpo_uni10_ctraction.setLabel([0,1])
    mpo_uni10_ctraction.setName('mpo_uni10_ctraction')
    mpo_uni10_6_list_prime[0]=mpo_uni10_ctraction

 ##printmpo_uni10_6_list_prime[0]

 ##############  Testing calcualtion of trace H^{2} ###########################
 mpo_uni10_br_2=copy.copy(mpo_boundy_list[3])
 mpo_uni10_br_2.setLabel([1])
 trH2=mpo_boundy_list[2]*mpo_uni10_6_list_prime[0]*mpo_uni10_br_2
 print trH2[0]
 ###############################################################################


 ##printmpo_uni10_6_list_prime[0].profile()
 ##printmpo_uni10_6_list_prime[0].printDiagram()
###################  intialize_unitary_list #################################### 
def intialize_unitary_list( L, L_lay,d, delta): 
 U_list=[]
 bdi_spin = uni10.Bond(uni10.BD_IN, d);
 bdo_spin = uni10.Bond(uni10.BD_OUT, d);
 Svd=[]
 matrix_random=uni10.Matrix(d*d, d*d)
 for i in xrange(d*d):
  for j in xrange(d*d):
   if(i==j):
    matrix_random[i*d*d+j]=1.0
   else:
    matrix_random[i*d*d+j]=0.0
 #matrix_random.setIdentity()
 for i in xrange(L):
  U_list.append([]) 
 for i in xrange(L):
  for j in xrange(len(L_lay)):
    U_uni10=uni10.UniTensor([bdi_spin, bdi_spin, bdo_spin,bdo_spin], "Unitary_uni10")
    matrix_random1=copy.copy(matrix_random)
    matrix_random1.randomize()
    matrix_random1=matrix_random+(delta)*matrix_random1
    Svd=matrix_random1.svd()
    U_uni10.putBlock( Svd[0]*Svd[2] )
    U_list[i].append(U_uni10)
 return U_list 


def unitary_list_transpose( U_list, L, L_lay):
  U_list_Trans=[]
  for i in xrange(L):
   U_list_Trans.append([]) 
  for i in xrange(L):
   for j in xrange(len(L_lay)):
    U_list_Trans[i].append(copy.copy(U_list[i][j]))
    U_list_Trans[i][j].transpose()
    U_list_Trans[i][j].setLabel([0,1,2,3])
  return U_list_Trans
    

def  make_mpo_U(U_list, L_position, L_lay, L, Letter):
 U=[copy.copy(U_list[L_position][i]) for i in xrange(len(L_lay))]
 bdo = uni10.Bond(uni10.BD_OUT, 1)
 bdi = uni10.Bond(uni10.BD_IN, 1)
 Uni_out_1=uni10.UniTensor([bdo], "Uni_out_1")
 Uni_in_1=uni10.UniTensor([bdo], "Uni_in_1")
 Uni_in_1.identity()
 Uni_out_1.identity()
 if (L_position==0):
  if (len(L_lay) == 1 ):
    U0=copy.copy(U[0])
    U0.combineBond([0,1])
    U0.combineBond([2,3])
    Uni_out_1.setLabel([4])
    U_result=U0*Uni_out_1
    U_result.permute([0,4,2],1)
    if Letter is 'up':
     U_result.setLabel([0,1,2])
     return U_result
    elif Letter is 'down':
     U_result.permute([2,4,0],1)
     U_result.setLabel([0,1,2])
     return U_result
  if len(L_lay) == 2 :
    U0=copy.copy(U[0])
    U1=copy.copy(U[1]) 
    U0.setLabel([0,1,2,3])
    U1.setLabel([3,4,5,6])
    U_result=U0*U1
    U_result.combineBond([0,1])
    U_result.combineBond([2,5,6])
    U_result.permute([0,4,2],1)
    U_result.setLabel([0,1,2])
    if Letter is 'up':
     return U_result
    elif Letter is 'down':
     U_result.permute([2,1,0],1)
     U_result.setLabel([0,1,2])
     return U_result
    
  if len(L_lay) == 3  :
    U0=copy.copy(U[0])
    U1=copy.copy(U[1])
    U2=copy.copy(U[2]) 
    U0.setLabel([0,1,2,3])
    U1.setLabel([3,4,5,6])
    U2.setLabel([2,5,7,8])
    U_result=((U0*U1)*(U2))
    U_result.combineBond([0,1])
    U_result.combineBond([7,8])
    U_result.combineBond([4,6])
    U_result.permute([0,4,7],1)
    U_result.setLabel([0,1,2])
    if Letter is 'up':
     return U_result
    elif Letter is 'down':
     U_result.permute([2,1,0],1)
     U_result.setLabel([0,1,2])
     return U_result
    
  if len(L_lay) == 4  :
    U0=copy.copy(U[0])
    U1=copy.copy(U[1])
    U2=copy.copy(U[2])
    U3=copy.copy(U[3]) 
    U0.setLabel([0,1,2,3])
    U1.setLabel([3,4,5,6])
    U2.setLabel([2,5,7,8])
    U3.setLabel([8,9,10,11])
    U_result=((U0*U1)*(U2*U3))
    U_result.combineBond([0,1])
    U_result.combineBond([7,10,11])
    U_result.combineBond([4,6,9])
    U_result.permute([0,4,7],1)
    U_result.setLabel([0,1,2])
    if Letter is 'up':
     return U_result
    elif Letter is 'down':
     U_result.permute([2,1,0],1)
     U_result.setLabel([0,1,2])
     return U_result
    
  if len(L_lay) == 5  :
    U0=copy.copy(U[0])
    U1=copy.copy(U[1])
    U2=copy.copy(U[2])
    U3=copy.copy(U[3]) 
    U4=copy.copy(U[4]) 
    U0.setLabel([0,1,2,3])
    U1.setLabel([3,4,5,6])
    U2.setLabel([2,5,7,8])
    U3.setLabel([8,9,10,11])
    U4.setLabel([7,10,12,13])
    U_result=(((U0*U1)*(U2*U3))*U4)
    U_result.combineBond([0,1])
    U_result.combineBond([12,13])
    U_result.combineBond([4,6,9,11])
    U_result.permute([0,4,12],1)
    U_result.setLabel([0,1,2])
    if Letter is 'up':
     return U_result
    elif Letter is 'down':
     U_result.permute([2,1,0],1)
     U_result.setLabel([0,1,2])
     return U_result
    
 if L_position==((L/2) - 1):
  if len(L_lay) == 1:
    U0=copy.copy(U[0])
    U0.combineBond([0,1])
    U0.combineBond([2,3])
    Uni_out_1.setLabel([4])
    U_result=U0*Uni_out_1
    U_result.permute([0,4,2],2)
    if (Letter is 'up'):
     U_result.setLabel([0,1,2])
     return U_result
    elif(Letter is 'down'):
     U_result.permute([2,4,0],2)
     U_result.setLabel([0,1,2])
     return U_result
     
  if len(L_lay) == 2  :
      U0=copy.copy(U[0])
      U0.combineBond([0,1])
      U0.permute([0,2,3],2)
      U0.setLabel([0,1,2])
      U_result=copy.copy(U0)
      if (Letter is 'up'):
       return U_result
      elif Letter is 'down':
       U_result.permute([2,1,0],2)
       U_result.setLabel([0,1,2])
       return U_result
     
  if len(L_lay) == 3  :
      U0=copy.copy(U[0])
      U0.combineBond([0,1])
      U2=copy.copy(U[2])
      U2.setLabel([4,3,5,6])
      U2.combineBond([5,6])
      U_result=U0*U2
      U_result.combineBond([2,4])
      U_result.permute([0,2,5],2)
      U_result.setLabel([0,1,2])
      ##printU_result.printDiagram()
      if Letter is 'up':
       return U_result
      elif Letter is 'down':
       U_result.permute([2,1,0],2)
       U_result.setLabel([0,1,2])
       return U_result
  
  if len(L_lay) is 4 :
      U0=copy.copy(U[0])
      U0.combineBond([0,1])
      U2=copy.copy(U[2])
      U2.setLabel([4,3,5,6])
      U_result=U0*U2
      U_result.combineBond([2,4,5])
      U_result.permute([0,2,6],2)
      U_result.setLabel([0,1,2])
      if Letter is 'up':
       return U_result
      elif Letter is 'down':
       U_result.permute([2,1,0],2)
       U_result.setLabel([0,1,2])
       return U_result
 
 
  if len(L_lay) is 5 :
      U0=copy.copy(U[0])
      U0.combineBond([0,1])
      U2=copy.copy(U[2])
      U2.setLabel([4,3,5,6])
      U4=copy.copy(U[4])
      U4.setLabel([7,6,8,9])
      U4.combineBond([8,9])
      U_result=((U0*U2)*(U4))
      U_result.combineBond([2,4,5,7])
      U_result.permute([0,2,8],2)
      U_result.setLabel([0,1,2])
      if Letter is 'up':
       return U_result       
      elif Letter is 'down':
       U_result.permute([2,1,0],2)
       U_result.setLabel([0,1,2])
       return U_result
      
 if  (L_position is not ((L/2) - 1)) and (L_position is not 0)  : 
  if len(L_lay) is 1  :
    U0=copy.copy(U[0])
    U0.combineBond([0,1])
    U0.combineBond([2,3])
    Uni_in_1.setLabel([5])
    Uni_out_1.setLabel([4])
    U_result=U0*Uni_in_1*Uni_out_1
    U_result.permute([5,0,4,2],2)
    U_result.setLabel([0,1,2,3])
    if Letter is 'up':
     return U_result
    elif Letter is 'down':
     U_result.permute([0,3,2,1],2)
     U_result.setLabel([0,1,2,3])
     return U_result


  if len(L_lay) is 2  :
    U0=copy.copy(U[0])
    U1=copy.copy(U[1]) 
    U0.setLabel([0,1,2,3])
    U1.setLabel([3,4,5,6])
    U_result=uni10.contract(U0,U1)
    U_result.combineBond([0,1])
    U_result.combineBond([5,6])
    U_result.permute([2,0,4,5],2)
    U_result.setLabel([0,1,2,3])
    #printU_result.printDiagram()
    if (Letter is 'up'):
     return U_result
    elif(Letter is 'down'):
     U_result.permute([0,3,2,1],2)
     U_result.setLabel([0,1,2,3])
     return U_result
    
  if len(L_lay) == 3 :
    U0=copy.copy(U[0])
    U1=copy.copy(U[1])
    U2=copy.copy(U[2]) 
    U0.setLabel([0,1,2,3])
    U1.setLabel([3,4,5,6])
    U2.setLabel([7,5,8,9])
    U_result=((U0*U1)*(U2))
    U_result.combineBond([0,1])
    U_result.combineBond([8,9])
    U_result.combineBond([4,6])
    U_result.combineBond([2,7])
    U_result.permute([2,0,4,8],2)
    U_result.setLabel([0,1,2,3])
    #printU_result.printDiagram()
    if Letter is 'up':
     return U_result
    elif Letter is 'down':
     U_result.permute([0,3,2,1],2)
     U_result.setLabel([0,1,2,3])
     return U_result
 
  if len(L_lay) == 4 :
    U0=copy.copy(U[0])
    U1=copy.copy(U[1])
    U2=copy.copy(U[2])
    U3=copy.copy(U[3])
    U0.setLabel([0,1,2,3])
    U1.setLabel([3,4,5,6])
    U2.setLabel([7,5,8,9])
    U3.setLabel([9,10,11,12])
    U_result=((U0*U1)*(U2*U3))
    U_result.combineBond([0,1])
    U_result.combineBond([11,12])
    U_result.combineBond([2,7,8])
    U_result.combineBond([4,6,10])
    U_result.permute([2,0,4,11],2)
    U_result.setLabel([0,1,2,3])
    if Letter is 'up':
     return U_result
    elif Letter is 'down':
     U_result.permute([0,3,2,1],2)
     U_result.setLabel([0,1,2,3])
     return U_result

  if (len(L_lay) == 5  ):
    U0=copy.copy(U[0])
    U1=copy.copy(U[1])
    U2=copy.copy(U[2])
    U3=copy.copy(U[3]) 
    U4=copy.copy(U[4]) 
    U0.setLabel([0,1,2,3])
    U1.setLabel([3,4,5,6])
    U2.setLabel([7,5,8,9])
    U3.setLabel([9,10,11,12])
    U4.setLabel([13,11,14,15])
    U_result=(((U0*U1)*(U2*U3))*(U4))
    U_result.combineBond([0,1])
    U_result.combineBond([14,15])
    U_result.combineBond([2,7,8,13])
    U_result.combineBond([4,6,10,12])
    U_result.permute([2,0,4,14],2)
    U_result.setLabel([0,1,2,3])
    #printU_result.printDiagram()
    if Letter is 'up':
     return U_result
    elif Letter is 'down':
     U_result.permute([0,3,2,1],2)
     U_result.setLabel([0,1,2,3])
     return U_result



def make_mpo_U_Label( L_position,L_lay, L, Letter):
 if (L_position==0):
  Bond_IN=2
  if (len(L_lay) == 1 ):
    if Letter=='up':
     per_labels=[0,1,4,2,3]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[2,3,4,0,1]
     return per_labels, Bond_IN
     
  if (len(L_lay) == 2 ):     
    if Letter=='up':
     per_labels=[0,1,4,2,5,6]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[2,5,6,4, 0,1]
     return per_labels, Bond_IN

  if (len(L_lay) == 3 ):
    if Letter=='up':
     per_labels=[0,1,4,6,7,8 ]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[4,6,7,8,0,1 ]
     return per_labels, Bond_IN

  if (len(L_lay) == 4 ):     
    if Letter=='up':
     per_labels=[0,1,4,6,9,7,10,11]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[7,10,11,4,6,9,0,1]
     return per_labels, Bond_IN

  if (len(L_lay) == 5 ):
    if Letter=='up':
     per_labels=[0,1,4,6,9,11,12,13]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[12,13,4,6,9,11,0,1]
     return per_labels, Bond_IN

 if L_position==((L/2)-1):
 
  if (len(L_lay) == 1 ):    
    Bond_IN=3
    if Letter=='up':
     per_labels=[0,1,4,2,3]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[2,3,4, 0,1]
     return per_labels, Bond_IN

  if (len(L_lay) == 2 ):
    Bond_IN=3  
    if Letter=='up':
     per_labels=[0,1,2,3]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[3,2,0,1]
     return per_labels, Bond_IN

  if (len(L_lay) == 3 ):
    Bond_IN=4    
    if Letter=='up':
     per_labels=[0,1,2,4,5,6]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[5,6,2,4,0,1]
     return per_labels, Bond_IN

  if (len(L_lay) == 4 ):
    Bond_IN=5    
    if Letter=='up':
     per_labels=[0,1,2,4,5,6]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[6,2,4,5,0,1]
     return per_labels, Bond_IN

  if (len(L_lay) == 5 ):
    Bond_IN=6    
    if Letter=='up':
     per_labels=[0,1,2,4,5,7,8,9]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[8,9,2,4,5,7,0,1]
     return per_labels, Bond_IN
 if  (L_position is not ((L/2) - 1)) and (L_position is not 0)  : 
  if (len(L_lay) == 1 ):
    Bond_IN=3    
    if Letter=='up':
     per_labels=[5,0,1,4,2,3]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[5,2,3,4,0,1]
     return per_labels, Bond_IN

  if (len(L_lay) == 2 ):
    Bond_IN=3    
    if Letter=='up':
     per_labels=[2,0,1,4,5,6]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[2,5,6,4,0,1]
     return per_labels, Bond_IN


  if (len(L_lay) == 3 ):
    Bond_IN=4    
    if Letter=='up':
     per_labels=[2,7,0,1,4,6,8,9]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[2,7,8,9,4,6,0,1]
     return per_labels, Bond_IN

  if (len(L_lay) == 4 ):
    Bond_IN=5    
    if Letter=='up':
     per_labels=[2,7,8,0,1,4,6,10,11,12]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[2,7,8,11,12,4,6,10,0,1]
     return per_labels, Bond_IN

  if (len(L_lay) == 5 ):
    Bond_IN=6    
    if Letter=='up':
     per_labels=[2,7,8,13,0,1,4,6,10,12,14,15]
     return per_labels, Bond_IN
    elif Letter=='down':
     per_labels=[2,7,8,13,14,15,4,6,10,12,0,1]
     return per_labels, Bond_IN


def make_mpo_U_list(U_list, L_lay, L, Letter):
 mpo_U_list=[]
 for i in xrange(L/2):
  mpo_U_list.append(make_mpo_U(U_list,i, L_lay, L, Letter))
 
 return mpo_U_list




