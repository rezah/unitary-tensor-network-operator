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
     



def Reduction_lastbond(U_up_1, n, Letter):
 Mat=U_up_1.getBlock()
 dim0=U_up_1.bond(0).dim()
 dim1=U_up_1.bond(1).dim()
 dim2=U_up_1.bond(2).dim()
 dim3=U_up_1.bond(3).dim()
 if Letter is 'up':
  Mat1=uni10.Matrix(dim0*dim1, dim2)
  Mat1.set_zero()
  for i in xrange(dim0):
   for j in xrange(dim1):
    for m in xrange(dim2):
     Mat1[i*dim1*dim2+j*dim2+m]=U_up_1[i*dim1*dim2*dim3+j*dim2*dim3+m*dim3+n]
  bond_list=list(U_up_1.bond())
  #print list(bond_list)
  bond_list=list(bond_list)
  bond_list.pop()
  #print bond_list
  Uni_tensor=uni10.UniTensor(bond_list,U_up_1.getName())
  Uni_tensor.putBlock(Mat1)
  #print Uni_tensor.printDiagram(), n
  return Uni_tensor
 elif Letter is 'down':
  Mat1=uni10.Matrix(dim0, dim2*dim3)
  Mat1.set_zero()
  for i in xrange(dim0):
   for j in xrange(dim3):
    for m in xrange(dim2):
     Mat1[i*dim3*dim2+m*dim3+j]=U_up_1[i*dim1*dim2*dim3+n*dim2*dim3+m*dim3+j]
  bond_list=list(U_up_1.bond())
  #print list(bond_list)
  bond_list=list(bond_list)
  A=bond_list.pop()
  B=bond_list.pop()
  C=bond_list.pop()
  bond_list.append(B)
  bond_list.append(A)
  #print bond_list
  Uni_tensor=uni10.UniTensor(bond_list,U_up_1.getName())
  Uni_tensor.putBlock(Mat1)
  #print Uni_tensor.printDiagram(), n
  return Uni_tensor
 
 
 
 
 
  
  
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
  print 'L_position=', L_position
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
  EnvLeft_net.putTensor('mpo_up',mpo_up)
  EnvLeft_net.putTensor('mpo_down',mpo_down)

  dim0=U_up_1.bond(3).dim()
  for i in xrange(dim0):
   U_up_11=Reduction_lastbond(U_up_1, i, 'up')
   U_down_11=Reduction_lastbond(U_down_1, i, 'down')
   U_up_22=Reduction_lastbond(U_up_2, i, 'up')
   U_down_22=Reduction_lastbond(U_down_2, i, 'down')
   EnvLeft_net.putTensor('U_up_11',U_up_11)
   EnvLeft_net.putTensor('U_up_22',U_up_22)
   EnvLeft_net.putTensor('U_down_11',U_down_11)
   EnvLeft_net.putTensor('U_down_22',U_down_22)
   if i is 0:
    Result_uni10=EnvLeft_net.launch()
    Result_uni10.setLabel([0,1,2,3,4,5])
    Result_final=copy.copy(Result_uni10)
    #print EnvLeft_net
   elif i is not 0:
    Result_uni10=EnvLeft_net.launch()
    Result_uni10.setLabel([0,1,2,3,4,5])
    Result_final=Result_final+Result_uni10
   
   
   #print Result_uni10.printDiagram()
  #print EnvLeft_net.profile()
  #print EnvLeft_net
  #print Result_uni10[5]
  return Result_final
  
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
  dim0=U_up_1.bond(3).dim()
  EnvRight_net = uni10.Network("EnvRight.net")
  EnvRight_net.putTensor('Env_right_copy', Env_right_copy)
  EnvRight_net.putTensor('mpo_up',mpo_up)
  EnvRight_net.putTensor('mpo_down',mpo_down)

  for i in xrange(dim0):
   U_up_11=Reduction_lastbond(U_up_1, i, 'up')
   U_down_11=Reduction_lastbond(U_down_1, i, 'down')
   U_up_22=Reduction_lastbond(U_up_2, i, 'up')
   U_down_22=Reduction_lastbond(U_down_2, i, 'down')
   EnvRight_net.putTensor('U_up_11',U_up_11)
   EnvRight_net.putTensor('U_up_22',U_up_22)
   EnvRight_net.putTensor('U_down_11',U_down_11)
   EnvRight_net.putTensor('U_down_22',U_down_22)
   if i is 0:
    Result_uni10=EnvRight_net.launch()
    Result_uni10.setLabel([0,1,2,3,4,5])
    Result_final=copy.copy(Result_uni10)
   elif i is not 0:
    Result_uni10=EnvRight_net.launch()
    Result_uni10.setLabel([0,1,2,3,4,5])
    Result_final=Result_final+Result_uni10

  return Result_final
  
  
def put_uni10(U_up_2,Result_final, Result_uni10, n):
  dim0=U_up_2.bond(0).dim()
  dim1=U_up_2.bond(1).dim()
  dim2=U_up_2.bond(2).dim()
  dim3=U_up_2.bond(3).dim()
  matH=Result_uni10.getBlock()
  for i in xrange(dim0):
   for j in xrange(dim1):
    for m in xrange(dim2):
     Result_final[i*dim1*dim2*dim3+j*dim2*dim3+m*dim3+n]=matH[i*dim1*dim2+j*dim2+m]


 
  
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
  
  
  dim0=U_up_2.bond(0).dim()
  dim1=U_up_2.bond(1).dim()
  dim2=U_up_2.bond(2).dim()
  dim3=U_up_2.bond(3).dim()
  Result_final=uni10.Matrix(dim0*dim1,dim2*dim3)
  Result_final.set_zero()
  EnvUni_net = uni10.Network("EnvUni.net")
  EnvUni_net.putTensor('Env_Left_copy',Env_Left_copy)
  EnvUni_net.putTensor('Env_Right_copy',Env_Right_copy)
  EnvUni_net.putTensor('mpo_up',mpo_up)
  EnvUni_net.putTensor('mpo_down',mpo_down)

  for i in xrange(dim3):
   U_down_11=Reduction_lastbond(U_down_1, i, 'down')
   U_up_22=Reduction_lastbond(U_up_2, i, 'up')
   U_down_22=Reduction_lastbond(U_down_2, i, 'down')
   EnvUni_net.putTensor('U_up_22',U_up_22)
   EnvUni_net.putTensor('U_down_11',U_down_11)
   EnvUni_net.putTensor('U_down_22',U_down_22)
   Result_uni10=EnvUni_net.launch()
   put_uni10(U_up_2,Result_final, Result_uni10, i)
   
  Result_final_uni=uni10.UniTensor(U_up_2.bond())
  Result_final_uni.putBlock(Result_final)
  Result_final_uni.setLabel([0,1,2,3])
  return Result_final_uni
 
 
def Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay, L_lay_selected, L_position, Env_Uni_inner, Tech):
 U=[copy.copy(U_list[L_position][i]) for i in xrange(len(L_lay))]
 if Tech is 'Regular': 
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
    Environment_Uni_newlabel.setLabel(perl_label_up[3])
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
    for i in xrange(len(perl_label_up[3])):
     if i < Bond_IN[3] :
      bd_list.append(bdi)
     else:
      bd_list.append(bdo)
    Environment_Uni_newlabel=uni10.UniTensor(bd_list)
    Environment_Uni_newlabel.setLabel(perl_label_up[3])
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
    for i in xrange(len(perl_label_up[3])):
     if i < Bond_IN[3] :
      bd_list.append(bdi)
     else:
      bd_list.append(bdo)
    Environment_Uni_newlabel=uni10.UniTensor(bd_list)
    Environment_Uni_newlabel.setLabel(perl_label_up[3])
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
    for i in xrange(len(perl_label_up[3])):
     if i < Bond_IN[3] :
      bd_list.append(bdi)
     else:
      bd_list.append(bdo)
    Environment_Uni_newlabel=uni10.UniTensor(bd_list)
    Environment_Uni_newlabel.setLabel(perl_label_up[3])
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
    for i in xrange(len(perl_label_up[3])):
     if i < Bond_IN[3] :
      bd_list.append(bdi)
     else:
      bd_list.append(bdo)
    Environment_Uni_newlabel=uni10.UniTensor(bd_list)
    Environment_Uni_newlabel.setLabel(perl_label_up[3])
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
  if  0 < L_position < (len(U_list)-1) :
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




 elif Tech is 'MERA':
  if L_position is 0:
   if len(L_lay) is 6:
    print 'hi'
    bdi=uni10.Bond(uni10.BD_IN, 4)
    bdo=uni10.Bond(uni10.BD_OUT, 4)
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
    U1.setLabel([2,4,5,6])
    if L_lay_selected is 0:
     U_result=Environment_Uni_newlabel*(U1)
     U_result.permute([0,1,2,3],2)
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 1:
     U_result=(Environment_Uni_newlabel*U0)
     U_result.permute([2,4,5,6],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
   if len(L_lay) is 8:
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
   if len(L_lay) is 6:
    bdi=uni10.Bond(uni10.BD_IN, 4)
    bdo=uni10.Bond(uni10.BD_OUT, 4)
    bd_list=[]
    for i in xrange(len(perl_label_up[3])):
     if i < Bond_IN[3]:
      bd_list.append(bdi)
     else:
      bd_list.append(bdo)
    Environment_Uni_newlabel=uni10.UniTensor(bd_list)
    Environment_Uni_newlabel.setLabel(perl_label_up[3])
    Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
    U0=copy.copy(U[0])
    U1=copy.copy(U[1])
    U2=copy.copy(U[2])
    U0.setLabel([0,1,2,3])
    U1.setLabel([4,2,5,6])
    U2.setLabel([6,3,7,8])


   if len(L_lay) is 8:
    bdi=uni10.Bond(uni10.BD_IN, 2)
    bdo=uni10.Bond(uni10.BD_OUT, 2)
    bd_list=[]
    for i in xrange(len(perl_label_up[3])):
     if i < Bond_IN[3] :
      bd_list.append(bdi)
     else:
      bd_list.append(bdo)
    Environment_Uni_newlabel=uni10.UniTensor(bd_list)
    Environment_Uni_newlabel.setLabel(perl_label_up[3])
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


  if  0 < L_position < (len(U_list)-1) :
   if len(L_lay) is 6:
    bdi=uni10.Bond(uni10.BD_IN, 2)
    bdo=uni10.Bond(uni10.BD_OUT, 2)
    bdi1=uni10.Bond(uni10.BD_IN, 4)
    bdo1=uni10.Bond(uni10.BD_OUT, 4)


    bd_list=[]
    for i in xrange(len(perl_label_up[1])):
     if i < Bond_IN[1] :
      bd_list.append(bdi)
     else:
      bd_list.append(bdo)
    bd_list.append(bdi1)
    bd_list.append(bdi1)
    bd_list.append(bdi)
    bd_list.append(bdi)
    bd_list.append(bdi)
    bd_list.append(bdi)
    bd_list.append(bdo1)
    bd_list.append(bdo1)
    bd_list.append(bdo1)
    bd_list.append(bdo1)
    
    
    
      
    
    Environment_Uni_newlabel=uni10.UniTensor(bd_list)
    Environment_Uni_newlabel.setLabel(perl_label_up[1])
    Environment_Uni_newlabel.putBlock(Environment_Uni[L_position].getBlock())
    U0=copy.copy(U[0])
    U1=copy.copy(U[1])
    U2=copy.copy(U[2])
    U3=copy.copy(U[3]) 
    U4=copy.copy(U[4])
    U5=copy.copy(U[5]) 
    U0.setLabel([0,1,4,8])
    U1.setLabel([2,3,9,10])
    U2.setLabel([5,4,6,7])
    U2.combineBond([6,7])
    U3.setLabel([8,9,11,12])
    U3.combineBond([11,12])
    U4.setLabel([6,11,13,14])
    U5.setLabel([14,15,16,17])
    if L_lay_selected is 0:
     U_result=(Environment_Uni_newlabel*U5)*((U1*U2*U3)*(U4))
     U_result.permute([0,1,4,8],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 1:
     U_result=(Environment_Uni_newlabel*U5)*((U0*U2*U3)*(U4))
     U_result.permute([2,3,9,10],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 2:
     U_result=(Environment_Uni_newlabel*U5)*((U1*U0*U3)*(U4))
     U_result.permute([5,4,6],2)
     Uni_h=uni10.UniTensor([bdi,bdi,bdo,bdo])
     Uni_h.putBlock(U_result.getBlock())
     Uni_h.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=Uni_h
    if L_lay_selected is 3:
     U_result=(Environment_Uni_newlabel*U5)*((U1*U2*U0)*(U4))
     U_result.permute([8,9,11],2)
     Uni_h=uni10.UniTensor([bdi,bdi,bdo,bdo])
     Uni_h.putBlock(U_result.getBlock())
     Uni_h.setLabel([0,1,2,3])     
     Env_Uni_inner[L_position][L_lay_selected]=Uni_h
    if L_lay_selected is 4:
     U_result=(Environment_Uni_newlabel*U5)*((U1*U3)*(U0*U2))
     U_result.permute([6,11,13,14],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 5:
     U_result=((Environment_Uni_newlabel)*(U0*U1*U2*U3))*(U4)
     U_result.permute([14,15,16,17],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result


   if (len(L_lay) is 8) and (L_position % 2 is 1) :
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
    U0.setLabel([0,1,2,3])
    U1.setLabel([4,5,6,7])
    U2.setLabel([8,9,10,11])
    if L_lay_selected is 0:
     U_result=Environment_Uni_newlabel*(U1*U2)
     U_result.permute([0,1,2,3],2)
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 1:
     U_result=Environment_Uni_newlabel*(U0*U2)
     U_result.permute([4,5,6,7],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 2:
     U_result=Environment_Uni_newlabel*(U0*U1)
     U_result.permute([8,9,10,11],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result

   if (len(L_lay) is 8) and (L_position % 2 is 0) :
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
    U1=copy.copy(U[1])
    U2=copy.copy(U[2])
    U3=copy.copy(U[3]) 
    U4=copy.copy(U[4]) 
    U5=copy.copy(U[5])
    U6=copy.copy(U[6])
    U7=copy.copy(U[7])
    U0.setLabel([0,1,2,3])
    U1.setLabel([4,2,5,6])
    U2.setLabel([3,7,8,9])
    U3.setLabel([6,8,10,11])
    U4.setLabel([12,13,14,15])
    U5.setLabel([5,14,17,18])
    U6.setLabel([15,9,19,16])
    U7.setLabel([18,19,20,21])
    if L_lay_selected is 0:
     U_result=(Environment_Uni_newlabel*(U4*U5*U6*U7))*(U1*U2*U3)
     U_result.permute([0,1,2,3],2)
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 1:
     U_result=(Environment_Uni_newlabel*(U4*U5*U6*U7))*(U0*U2*U3)
     U_result.permute([4,2,5,6],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 2:
     U_result=(Environment_Uni_newlabel*(U4*U5*U6*U7))*(U0*U1*U3)
     U_result.permute([3,7,8,9],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 3:
     U_result=(Environment_Uni_newlabel*(U4*U5*U6*U7))*(U0*U1*U2)
     U_result.permute([6,8,10,11],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 4:
     U_result=(Environment_Uni_newlabel*(U0*U1*U2*U3))*(U5*U6*U7)
     U_result.permute([12,13,14,15],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 5:
     U_result=(Environment_Uni_newlabel*(U0*U1*U2*U3))*(U4*U6*U7)
     U_result.permute([5,14,17,18],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 6:
     U_result=(Environment_Uni_newlabel*(U0*U1*U2*U3))*(U5*U4*U7)
     U_result.permute([15,9,19,16],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result
    if L_lay_selected is 7:
     U_result=(Environment_Uni_newlabel*(U0*U1*U2*U3))*(U4*U5*U6)
     U_result.permute([18,19,20,21],2)
     U_result.setLabel([0,1,2,3])
     Env_Uni_inner[L_position][L_lay_selected]=U_result

   
   
   
   
