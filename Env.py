import pyUni10 as uni10
#import sys
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import random 
import copy

def Iden_uni10(d):
 bdi=uni10.Bond(uni10.BD_IN, d)
 bdo=uni10.Bond(uni10.BD_OUT, d)
 Iden_Matrix=uni10.Matrix(d*d,d)
 for i in xrange(d):
   for j in xrange(d):
    for m in xrange(d):
     if i is j and j is m: 
      Iden_Matrix[i*d+j,m]=1
     else:
      Iden_Matrix[i*d+j,m]=0
 Iden_uni10=uni10.UniTensor([bdi,bdi,bdo], 'Iden_uni10')
 Iden_uni10.putBlock(Iden_Matrix)   
 return Iden_uni10
     


def Env_left ( mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left   ):
 
 if L_position is 0:
  U_up_1=copy.copy(mpo_U_list_up[L_position])
  U_down_1=copy.copy(mpo_U_list_down[L_position])
  U_up_2=copy.copy(mpo_U_list_up[L_position])
  U_down_2=copy.copy(mpo_U_list_down[L_position])
  mpo_up=copy.copy(mpo_list2[L_position])
  mpo_down=copy.copy(mpo_list2[L_position])
  mpo_boundy_up=copy.copy(mpo_boundy_list[L_position])
  mpo_boundy_down=copy.copy(mpo_boundy_list[L_position])
  Iden_up=Iden_uni10( U_up_1.bond(2).dim() )
  Iden_down=Iden_uni10(U_up_1.bond(2).dim())
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
  ##########  Contraction  ##############
  A_uni10=(((U_up_1*mpo_up)*(mpo_boundy_up))*(U_down_1))*Iden_up
  B_uni10=(((U_up_2*mpo_down)*(mpo_boundy_down))*(U_down_2))*Iden_down
  Environment_Left=A_uni10*B_uni10
  Environment_Left.permute([0,1,2,3,4,5],0)
  Environment_Left.setName('Environment_Left')
  #print Environment_Left.printDiagram()
  #print Environment_Left
  return Environment_Left
 else: 
  Env_left_copy=copy.copy(Environment_Left[L_position-1]) 
  U_up_1=copy.copy(mpo_U_list_up[L_position])
  U_down_1=copy.copy(mpo_U_list_down[L_position])
  U_up_2=copy.copy(mpo_U_list_up[L_position])
  U_down_2=copy.copy(mpo_U_list_down[L_position])
  mpo_up=copy.copy(mpo_list2[L_position])
  mpo_down=copy.copy(mpo_list2[L_position])
  Iden_up=Iden_uni10(U_up_1.bond(3).dim())
  Iden_down=Iden_uni10(U_up_1.bond(3).dim())
#  ##### assign label#####
#  U_up_1.setLabel([0,6,-1,8])
#  U_up_2.setLabel([3,11,-4,14])
#  U_down_1.setLabel([2,9,-3,7])
#  U_down_2.setLabel([5,13,-6,12])
#  mpo_up.setLabel([1,7,-2,6])
#  mpo_down.setLabel([4,12,-5,11])
#  Iden_up.setLabel([10,9,8])
#  Iden_down.setLabel([10,13,14])
#  ###########   Contraction  ####################
#  A_uni10=((Environment_Left[L_position-1]*(U_up_1*mpo_up))*U_down_1)*Iden_up
#  Result_uni10=((A_uni10*(U_up_2*mpo_down))*U_down_2)*Iden_down
#  print Result_uni10.printDiagram()
#  Result_uni10.permute([-1,-2,-3,-4,-5,-6],0)
#  Result_uni10.setLabel([0,1,2,3,4,5])
#  Result_uni10.setName('Environment_Left')
#  print Result_uni10.profile()
#  print Result_uni10[5]
  
  ############## Contraction using Network ###############
  EnvLeft_net = uni10.Network("EnvLeft.net")
  EnvLeft_net.putTensor('Env_left_copy',Env_left_copy)
  EnvLeft_net.putTensor('U_up_1',U_up_1)
  EnvLeft_net.putTensor('U_up_2',U_up_2)
  EnvLeft_net.putTensor('U_down_1',U_down_1)
  EnvLeft_net.putTensor('U_down_2',U_down_2)
  EnvLeft_net.putTensor('mpo_up',mpo_up)
  EnvLeft_net.putTensor('mpo_down',mpo_down)
  EnvLeft_net.putTensor('Iden_up',Iden_up)
  EnvLeft_net.putTensor('Iden_down',Iden_down)
  Result_uni10=EnvLeft_net.launch()
  Result_uni10.setLabel([0,1,2,3,4,5])
  #print Result_uni10.printDiagram()
  #print EnvLeft_net.profile()
  #print EnvLeft_net
  #print Result_uni10[5]
  return Result_uni10
  
def Env_right ( mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Right   ):
 
 if L_position is (len(mpo_list2)-1):
  U_up_1=copy.copy(mpo_U_list_up[L_position])
  U_down_1=copy.copy(mpo_U_list_down[L_position])
  U_up_2=copy.copy(mpo_U_list_up[L_position])
  U_down_2=copy.copy(mpo_U_list_down[L_position])
  mpo_up=copy.copy(mpo_list2[L_position])
  mpo_down=copy.copy(mpo_list2[L_position])
  mpo_boundy_up=copy.copy(mpo_boundy_list[1])
  mpo_boundy_down=copy.copy(mpo_boundy_list[1])
  Iden_up=Iden_uni10(U_up_1.bond(2).dim())
  Iden_down=Iden_uni10(U_up_1.bond(2).dim())
  ##### assign label#####
  U_up_1.setLabel([-1,0,-5])
  U_up_2.setLabel([-9,3,-7])
  U_down_1.setLabel([-4,2,-3])
  U_down_2.setLabel([-8,5,-11])
  mpo_up.setLabel([1,-3,-2,-1])
  mpo_down.setLabel([4,-11,-10,-9])
  mpo_boundy_up.setLabel([-2])
  mpo_boundy_down.setLabel([-10])
  Iden_up.setLabel([-6,-4,-5])
  Iden_down.setLabel([-6,-8,-7])
  ##########  Contraction  ##############
  A_uni10=(((U_up_1*mpo_up)*(mpo_boundy_up))*(U_down_1))*Iden_up
  B_uni10=(((U_up_2*mpo_down)*(mpo_boundy_down))*(U_down_2))*Iden_down
  Environment_Right=A_uni10*B_uni10
  Environment_Right.permute([0,1,2,3,4,5],6)
  Environment_Right.setName('Environment_Right')
  #print Environment_Right.printDiagram()
  #print Environment_Right
  return Environment_Right
 else: 
  Env_right_copy=copy.copy(Environment_Right[L_position+1]) 
  U_up_1=copy.copy(mpo_U_list_up[L_position])
  U_down_1=copy.copy(mpo_U_list_down[L_position])
  U_up_2=copy.copy(mpo_U_list_up[L_position])
  U_down_2=copy.copy(mpo_U_list_down[L_position])
  mpo_up=copy.copy(mpo_list2[L_position])
  mpo_down=copy.copy(mpo_list2[L_position])
  Iden_up=Iden_uni10(U_up_1.bond(3).dim())
  Iden_down=Iden_uni10(U_up_1.bond(3).dim())
  ############## Contraction using Network ###############
  EnvRight_net = uni10.Network("EnvRight.net")
  EnvRight_net.putTensor('Env_right_copy', Env_right_copy)
  EnvRight_net.putTensor('U_up_1',U_up_1)
  EnvRight_net.putTensor('U_up_2',U_up_2)
  EnvRight_net.putTensor('U_down_1',U_down_1)
  EnvRight_net.putTensor('U_down_2',U_down_2)
  EnvRight_net.putTensor('mpo_up',mpo_up)
  EnvRight_net.putTensor('mpo_down',mpo_down)
  EnvRight_net.putTensor('Iden_up',Iden_up)
  EnvRight_net.putTensor('Iden_down',Iden_down)
  Result_uni10=EnvRight_net.launch()
  Result_uni10.setLabel([0,1,2,3,4,5])
  #print Result_uni10.printDiagram()
  #print EnvRight_net.profile()
  #print EnvRight_net
  #print Result_uni10[5]
  return Result_uni10
  
  
  
  
def Environment_uni(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right):
 if L_position is 0:
  Environment_Right_copy=copy.copy(Environment_Right[L_position+1])
  U_down_1=copy.copy(mpo_U_list_down[0])
  U_up_2=copy.copy(mpo_U_list_up[0])
  U_down_2=copy.copy(mpo_U_list_down[0])
  mpo_up=copy.copy(mpo_list2[0])
  mpo_down=copy.copy(mpo_list2[0])
  mpo_boundy_up=copy.copy(mpo_boundy_list[0])
  mpo_boundy_down=copy.copy(mpo_boundy_list[0])
  Iden_up=Iden_uni10(U_up_2.bond(2).dim())
  Iden_down=Iden_uni10(U_up_2.bond(2).dim())
  ##### assign label#####
  U_up_2.setLabel([-9,3,-7])
  U_down_1.setLabel([-4,2,-3])
  U_down_2.setLabel([-8,5,-11])
  mpo_up.setLabel([-2,-3,1,-1])
  mpo_down.setLabel([-10,-11,4,-9])
  mpo_boundy_up.setLabel([-2])
  mpo_boundy_down.setLabel([-10])
  Iden_up.setLabel([-6,-4,-5])
  Iden_down.setLabel([-6,-8,-7])
  ##########  Contraction  ##############
  Environment_Uni=(((Environment_Right_copy*(mpo_boundy_up*mpo_up*U_down_1))*Iden_up)*(U_up_2*mpo_down*mpo_boundy_down))*(U_down_2*Iden_down)
  Environment_Uni.permute([-1,0,-5],1)
  Environment_Uni.setName('Environment_Uni')
  Environment_Uni.setLabel([0,1,2])
  print Environment_Uni.printDiagram()
  #print Environment_Left
  #return Environment_Left

 
 
 
 
 
 
 
  
  
