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
  
  
  
  
def Environment_uni_function(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right):
 if L_position is 0:
  Env_Right_copy=copy.copy(Environment_Right[L_position+1])
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
  Environment_Uni=(((Env_Right_copy*(mpo_boundy_up*mpo_up*U_down_1))*Iden_up)*(U_up_2*mpo_down*mpo_boundy_down))*(U_down_2*Iden_down)
  Environment_Uni.permute([-1,0,-5],1)
  Environment_Uni.setName('Environment_Uni')
  Environment_Uni.setLabel([0,1,2])
  #print Environment_Uni.printDiagram()
  #print Environment_Uni
  return Environment_Uni
 elif L_position is (len(mpo_list2)-1):
  Env_Left_copy=copy.copy(Environment_Left[L_position-1])
  U_down_1=copy.copy(mpo_U_list_down[L_position])
  U_up_2=copy.copy(mpo_U_list_up[L_position])
  U_down_2=copy.copy(mpo_U_list_down[L_position])
  mpo_up=copy.copy(mpo_list2[L_position])
  mpo_down=copy.copy(mpo_list2[L_position])
  mpo_boundy_up=copy.copy(mpo_boundy_list[1])
  mpo_boundy_down=copy.copy(mpo_boundy_list[1])
  Iden_up=Iden_uni10(U_up_2.bond(2).dim())
  Iden_down=Iden_uni10(U_up_2.bond(2).dim())
  ##### assign label#####
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
  Environment_Uni=(((((Env_Left_copy)*(mpo_up*mpo_boundy_up))*(U_down_1))*Iden_up)*(U_up_2*mpo_down*mpo_boundy_down))*(U_down_2*Iden_down)
  Environment_Uni.permute([-1,0,-5],2)
  Environment_Uni.setName('Environment_Uni')
  Environment_Uni.setLabel([0,1,2]) 
  #print Environment_Uni.printDiagram()
  #print Environment_Uni
  return Environment_Uni
 else:
  Env_Right_copy=copy.copy(Environment_Right[L_position+1])
  Env_Left_copy=copy.copy(Environment_Left[L_position-1])
  U_down_1=copy.copy(mpo_U_list_down[L_position])
  U_up_2=copy.copy(mpo_U_list_up[L_position])
  U_down_2=copy.copy(mpo_U_list_down[L_position])
  mpo_up=copy.copy(mpo_list2[L_position])
  mpo_down=copy.copy(mpo_list2[L_position])
  Iden_up=Iden_uni10(U_up_2.bond(3).dim())
  Iden_down=Iden_uni10(U_up_2.bond(3).dim())
  ############## Contraction using Network ###############
  EnvUni_net = uni10.Network("EnvUni.net")
  EnvUni_net.putTensor('Env_Left_copy',Env_Left_copy)
  EnvUni_net.putTensor('Env_Right_copy',Env_Right_copy)
  EnvUni_net.putTensor('U_up_2',U_up_2)
  EnvUni_net.putTensor('U_down_1',U_down_1)
  EnvUni_net.putTensor('U_down_2',U_down_2)
  EnvUni_net.putTensor('mpo_up',mpo_up)
  EnvUni_net.putTensor('mpo_down',mpo_down)
  EnvUni_net.putTensor('Iden_up',Iden_up)
  EnvUni_net.putTensor('Iden_down',Iden_down)
  Result_uni10=EnvUni_net.launch()
  Result_uni10.setLabel([0,1,2,3])
  return Result_uni10

 
 
 
 
def Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay, L_lay_selected, L_position, Env_Uni_inner):
 U=[copy.copy(U_list[L_position][i]) for i in xrange(len(L_lay))]
 if L_position is 0:
  if len(L_lay) is 1:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd1=uni10.Bond(uni10.BD_OUT, 1)
   Environment_Uni_newlabel=uni10.UniTensor([bdi, bdi, bd1, bdo, bdo])
   Environment_Uni_newlabel.setLabel(perl_label_up[0])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   #print Environment_Uni_newlabel.printDiagram()
   Uni1=uni10.UniTensor([bd1])
   Uni1.identity()
   Uni1.setLabel([4])
   Environment_Uni_newlabel=Environment_Uni_newlabel*Uni1
   Environment_Uni_newlabel.permute([0,1, 2, 3],2)
   #print Environment_Uni_newlabel.printDiagram()
   Env_Uni_inner[L_position][L_lay_selected]=Environment_Uni_newlabel
  if len(L_lay) is 2:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[L_position])):
    if i < Bond_IN[L_position] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[L_position])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   U0=copy.copy(U[0])
   U1=copy.copy(U[1]) 
   U0.setLabel([0,1,2,3])
   U1.setLabel([3,4,5,6])
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*U1
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 1:
    U_result=Environment_Uni_newlabel*U0
    U_result.permute([3,4,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
  if len(L_lay) is 3:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[L_position])):
    if i < Bond_IN[L_position] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[L_position])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   U0=copy.copy(U[0])
   U1=copy.copy(U[1])
   U2=copy.copy(U[2]) 
   U0.setLabel([0,1,2,3])
   U1.setLabel([3,4,5,6])
   U2.setLabel([2,5,7,8])
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*(U1*U2)
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 1:
    U_result=(Environment_Uni_newlabel*U0)*U2
    U_result.permute([3,4,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 2:
    U_result=Environment_Uni_newlabel*(U0*U1)
    U_result.permute([2,5,7,8],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
  if len(L_lay) is 4:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[L_position])):
    if i < Bond_IN[L_position] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[L_position])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   U0=copy.copy(U[0])
   U1=copy.copy(U[1])
   U2=copy.copy(U[2])
   U3=copy.copy(U[3]) 
   U0.setLabel([0,1,2,3])
   U1.setLabel([3,4,5,6])
   U2.setLabel([2,5,7,8])
   U3.setLabel([8,9,10,11])
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*(U1*U2*U3)
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 1:
    U_result=((Environment_Uni_newlabel*U0)*U2)*U3
    U_result.permute([3,4,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 2:
    U_result=(Environment_Uni_newlabel*(U0*U1))*U3
    U_result.permute([2,5,7,8],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 3:
    U_result=(Environment_Uni_newlabel*(U0*U1*U2))
    U_result.permute([8,9,10,11],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
  if len(L_lay) is 5:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[L_position])):
    if i < Bond_IN[L_position] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[L_position])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
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
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*(U1*U2*U3*U4)
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 1:
    U_result=((Environment_Uni_newlabel*U0)*U2)*(U3*U4)
    U_result.permute([3,4,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 2:
    U_result=(Environment_Uni_newlabel*(U0*U1))*(U3*U4)
    U_result.permute([2,5,7,8],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 3:
    U_result=(Environment_Uni_newlabel*(U0*U1*U2))*U4
    U_result.permute([8,9,10,11],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 4:
    U_result=(Environment_Uni_newlabel*(U0*U1*U2*U3))
    U_result.permute([7,10,12,13],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
 if L_position is (len(U_list)-1):
  if len(L_lay) is 1:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd1=uni10.Bond(uni10.BD_IN, 1)
   Environment_Uni_newlabel=uni10.UniTensor([bdi, bdi, bd1, bdo, bdo])
   Environment_Uni_newlabel.setLabel(perl_label_up[2])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   #print Environment_Uni_newlabel.printDiagram()
   Uni1=uni10.UniTensor([bd1])
   Uni1.identity()
   Uni1.setLabel([4])
   Environment_Uni_newlabel=Environment_Uni_newlabel*Uni1
   Environment_Uni_newlabel.permute([0,1, 2, 3],2)
   #print Environment_Uni_newlabel.printDiagram()
   Env_Uni_inner[L_position][L_lay_selected]=Environment_Uni_newlabel
  if len(L_lay) is 2:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[2])):
    if i < Bond_IN[2] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[2])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   U0=copy.copy(U[0]) 
   U0.setLabel([0,1,2,3])
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   else : print 'L_selected is beyond, there is no Unitary'
  if len(L_lay) is 3:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[2])):
    if i < Bond_IN[2] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[2])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   U0=copy.copy(U[0])
   U0.setLabel([0,1,2,3])
   U2=copy.copy(U[2])
   U2.setLabel([4,3,5,6])
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*U2
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   elif L_lay_selected is 2:
    U_result=Environment_Uni_newlabel*U0
    U_result.permute([4,3,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   else: print 'there is no Unitary with L_selected'
  if len(L_lay) is 4:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[2])):
    if i < Bond_IN[2] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[2])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   U0=copy.copy(U[0])
   U0.setLabel([0,1,2,3])
   U2=copy.copy(U[2])
   U2.setLabel([4,3,5,6])
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*U2
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   elif L_lay_selected is 2:
    U_result=Environment_Uni_newlabel*U0
    U_result.permute([4,3,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   else: print 'there is no Unitary with L_selected'   
  if len(L_lay) is 5:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[2])):
    if i < Bond_IN[2] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[2])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   U0=copy.copy(U[0])
   U2=copy.copy(U[2])
   U4=copy.copy(U[4])
   U0.setLabel([0,1,2,3])
   U2.setLabel([4,3,5,6])
   U4.setLabel([7,6,8,9])
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*(U2*U4)
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   elif L_lay_selected is 2:
    U_result=(Environment_Uni_newlabel*U0)*U4
    U_result.permute([4,3,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   elif L_lay_selected is 4:
    U_result=Environment_Uni_newlabel*(U0*U2)
    U_result.permute([7,6,8,9],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   else: print 'there is no Unitary with L_selected'   
 if  0<L_position < (len(U_list)-1) :
  if len(L_lay) is 1:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bdi_1=uni10.Bond(uni10.BD_IN, 1)
   bdo_1=uni10.Bond(uni10.BD_OUT, 1)
   Environment_Uni_newlabel=uni10.UniTensor([bdi_1, bdi, bdi, bdo_1, bdo, bdo])
   Environment_Uni_newlabel.setLabel(perl_label_up[1])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   #print Environment_Uni_newlabel.printDiagram()
   Uni_in=uni10.UniTensor([bdi_1])
   Uni_out=uni10.UniTensor([bdo_1])
   Uni_in.identity()
   Uni_out.identity()
   Uni_in.setLabel([5])
   Uni_out.setLabel([4])
   Environment_Uni_newlabel=Environment_Uni_newlabel*Uni_in*Uni_out
   Environment_Uni_newlabel.permute([0,1, 2, 3],2)
   #print Environment_Uni_newlabel.printDiagram()
   Env_Uni_inner[L_position][L_lay_selected]=Environment_Uni_newlabel
  if len(L_lay) is 2:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[1])):
    if i < Bond_IN[1] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[1])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   U0=copy.copy(U[0])
   U1=copy.copy(U[1]) 
   U0.setLabel([0,1,2,3])
   U1.setLabel([3,4,5,6])
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*U1
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 1:
    U_result=Environment_Uni_newlabel*U0
    U_result.permute([3,4,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
  if len(L_lay) is 3:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[1])):
    if i < Bond_IN[1]:
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[1])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   U0=copy.copy(U[0])
   U1=copy.copy(U[1])
   U2=copy.copy(U[2]) 
   U0.setLabel([0,1,2,3])
   U1.setLabel([3,4,5,6])
   U2.setLabel([7,5,8,9])
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*(U1*U2)
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 1:
    U_result=(Environment_Uni_newlabel*U0)*U2
    U_result.permute([3,4,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 2:
    U_result=Environment_Uni_newlabel*(U0*U1)
    U_result.permute([7,5,8,9],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
  if len(L_lay) is 4:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[1])):
    if i < Bond_IN[1] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[1])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
   U0=copy.copy(U[0])
   U1=copy.copy(U[1])
   U2=copy.copy(U[2])
   U3=copy.copy(U[3]) 
   U0.setLabel([0,1,2,3])
   U1.setLabel([3,4,5,6])
   U2.setLabel([7,5,8,9])
   U3.setLabel([9,10,11,12])
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*(U1*U2*U3)
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 1:
    U_result=((Environment_Uni_newlabel*U0)*U2)*U3
    U_result.permute([3,4,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 2:
    U_result=(Environment_Uni_newlabel*(U0*U1))*U3
    U_result.permute([7,5,8,9],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 3:
    U_result=(Environment_Uni_newlabel*(U0*U1*U2))
    U_result.permute([9,10,11,12],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
  if len(L_lay) is 5:
   bdi=uni10.Bond(uni10.BD_IN, 2)
   bdo=uni10.Bond(uni10.BD_OUT, 2)
   bd_list=[]
   for i in xrange(len(perl_label_up[1])):
    if i < Bond_IN[1] :
     bd_list.append(bdi)
    else:
     bd_list.append(bdo)
   Environment_Uni_newlabel=uni10.UniTensor(bd_list)
   Environment_Uni_newlabel.setLabel(perl_label_up[1])
   Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
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
   if L_lay_selected is 0:
    U_result=Environment_Uni_newlabel*(U1*U2*U3*U4)
    U_result.permute([0,1,2,3],2)
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 1:
    U_result=((Environment_Uni_newlabel*U0)*U2)*(U3*U4)
    U_result.permute([3,4,5,6],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 2:
    U_result=(Environment_Uni_newlabel*(U0*U1))*(U3*U4)
    U_result.permute([7,5,8,9],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 3:
    U_result=(Environment_Uni_newlabel*(U0*U1*U2))*U4
    U_result.permute([9,10,11,12],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result
   if L_lay_selected is 4:
    U_result=(Environment_Uni_newlabel*(U0*U1*U2*U3))
    U_result.permute([13,11,14,15],2)
    U_result.setLabel([0,1,2,3])
    Env_Uni_inner[L_position][L_lay_selected]=U_result

