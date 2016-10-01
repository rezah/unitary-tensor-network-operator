import pyUni10 as uni10
#import sys
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import random 
import copy

def matS0():
  return uni10.Matrix(2, 2, [1.0, 0.0, 0.0, 1.0])

def matSx():
  return uni10.Matrix(2, 2, [0.0, 0.50, 0.50, 0.0])

def matSy():
  return uni10.Matrix(2, 2, [0.0, -0.50, 0.50, 0.0])

def matSz():
  return uni10.Matrix(2, 2, [0.50, 0.0, 0.0, -0.50])


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
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=sx[s*d+sp]
        if(mu==4 and mup==2  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=-sy[s*d+sp]
        if(mu==4 and mup==3  ):
           mpo_matrix[mu*d*d*chi+s*d*chi+mup*d+sp]=sz[s*d+sp]
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
   
 #print mpo_uni10_1_list[1].printDiagram()
 return mpo_uni10_1_list, mpo_uni10_4_list, mpo_uni10_bl , mpo_uni10_br, mpo_uni10_bl_2, mpo_uni10_br_2 






def contraction_MPO_trH2(mpo_uni10_1_list, mpo_uni10_4_list,  mpo_uni10_bl, mpo_uni10_br, mpo_uni10_bl_2, mpo_uni10_br_2):
 mpo_uni10_2_list=[]
 mpo_uni10_3_list=[]
 mpo_uni10_5_list=[]
 mpo_uni10_6_list=[]
 
 mpo_uni10_4_listN=[copy.copy(mpo_uni10_4_list[i]) for i in xrange(len(mpo_uni10_4_list)) ]
 mpo_uni10_1_listN=[copy.copy(mpo_uni10_1_list[i])  for i in xrange(len(mpo_uni10_1_list)) ]
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

 #print mpo_uni10_list[0], mpo_uni10_list[1]
 #print mpo_uni10_2_list[0], mpo_uni10_2_list[1]
 #print mpo_uni10_3_list[0], mpo_uni10_3_list[1]

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
 mpo_uni10_br_2.setLabel([1])
 trH2=mpo_uni10_bl_2*mpo_uni10_3_list_prime[0]*mpo_uni10_br_2
 print trH2[0]
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
  #print mpo_uni10_5_list[0], mpo_uni10_6_list[0] 
  



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

 #print mpo_uni10_6_list_prime[0]

 ##############  Testing calcualtion of trace H^{2} ###########################
 mpo_uni10_br_2.setLabel([1])
 trH2=mpo_uni10_bl_2*mpo_uni10_6_list_prime[0]*mpo_uni10_br_2
 print trH2[0]
 ###############################################################################


 #print mpo_uni10_6_list_prime[0].profile()
 #print mpo_uni10_6_list_prime[0].printDiagram()

