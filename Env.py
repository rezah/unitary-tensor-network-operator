import pyUni10 as uni10
#import sys
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import random 
import copy

def Iden_uni10(d):
 
 bdi=uni10.Bond(uni10.BD_IN, d*d)
 bdo=uni10.Bond(uni10.BD_OUT, d*d)
 Iden_Matrix=uni10.Matrix(d*d*d*d,d*d)
 for i in xrange(d*d):
   for j in xrange(d*d):
    for m in xrange(d*d):
     if i is j and j is m: 
      Iden_Matrix[i*d*d+j,m]=1
     else:
      Iden_Matrix[i*d*d+j,m]=0
 Iden_uni10=uni10.UniTensor([bdi,bdi,bdo], 'Iden_uni10')
 Iden_uni10.putBlock(Iden_Matrix)   
 return Iden_uni10
     


def Env_left ( mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d   ):
 
 if L_position is 0:
  U_up_1=copy.copy(mpo_U_list_up[0])
  U_down_1=copy.copy(mpo_U_list_down[0])
  U_up_2=copy.copy(mpo_U_list_up[0])
  U_down_2=copy.copy(mpo_U_list_down[0])
  mpo_up=copy.copy(mpo_list2[0])
  mpo_down=copy.copy(mpo_list2[0])
  mpo_boundy_up=copy.copy(mpo_boundy_list[0])
  mpo_boundy_down=copy.copy(mpo_boundy_list[0])
  Iden_up=Iden_uni10(d)
  Iden_down=Iden_uni10(d)
  
  ##### assign label#####
  U_up_1.setLabel([-1,0,-5])
  U_up_2.setLabel([-9,3,-7])
  U_down_1.setLabel([-4,2,-3])
  U_down_2.setLabel([-8,5,-11])
  mpo_up.setLabel([-2,-3,1,-1])
  mpo_down.setLabel([-10,-11,4,-9])
  mpo_boundy_up.setLabel([-2])
  mpo_boundy_down.setLabel([-10])
  Iden_up.setLabel([-6,-4,-5])
  Iden_down.setLabel([-6,-8,-7])
  ########################
  A_uni10=(((U_up_1*mpo_up)*(mpo_boundy_up))*(U_down_1))*Iden_up
  B_uni10=(((U_up_2*mpo_down)*(mpo_boundy_down))*(U_down_2))*Iden_down
  Environment_Left=A_uni10*B_uni10
  Environment_Left.permute([0,1,2,3,4,5],0)
  Environment_Left.setName('Environment_Left')
  #print Environment_Left.printDiagram()
  #print Environment_Left
  return Environment_Left
  
  
  
  
  
  
  
