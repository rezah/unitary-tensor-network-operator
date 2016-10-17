import pyUni10 as uni10
#import sys
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
import random 
import copy
import env
import mpo 
####################################      Energy           #####################################################
def Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN):
 
 Update_Unitary_change(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L)
 
 Environment_Uni[L_position]=env.Environment_uni_function(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right)
 
 env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner )
 E=(Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0]
 return E

#####################################     Update mpo_U_up/down & Environment_Left  ################################
def Update_Unitary_change(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L):
 mpo_U_up=mpo.make_mpo_U(U_list, L_position, L_lay, L, 'up')
 mpo_U_down=mpo.make_mpo_U(U_list, L_position, L_lay, L, 'down')
 mpo_U_list_up[L_position]=mpo_U_up
 mpo_U_list_down[L_position]=mpo_U_down
# Environment_Left[L_position]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
################################   Copy U_list       ############################################################
def copy_U_list_function(U_list,L,L_lay):
 U_list_copy=[]
 for i in xrange(L/2):
  U_list_copy.append([])
 for i in xrange(L/2):
  for j in xrange(len(L_lay)):
    U_list_copy[i].append([None])
 for i in xrange(L/2):
   for j in xrange(len(L_lay)):
    U_list_copy[i][j]=copy.copy(U_list[i][j])
 return U_list_copy
#############################################################################################

###################################Optimze  ###########################################################
def optimize_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin):
 
 #U_list_copy=copy_U_list_function(U_list,L,L_lay)
 
# for i in xrange(len(L_lay)):
#  L_lay_selected=i
#  optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method)

 if L_position is not L/2-1:
  for i in xrange(len(L_lay)):
   L_lay_selected=i
   optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method,Max_SVD_iteratoin, Max_Steepest_iteratoin)
 elif L_position is L/2 -1:
  for i in xrange(0,len(L_lay),2):
   L_lay_selected=i
   optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin)


 Energy_f=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN)

 #print 'E_f', Energy_f
 

################################ Optimize inner function #########################################
def optimize_inner_function(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method,Max_SVD_iteratoin, Max_Steepest_iteratoin):
 
 Energy_s=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN)
 #print 'E_s=', Energy_s
 
 if Method is 'SVD':
  U_update=copy.copy(U_list[L_position][L_lay_selected]).getBlock()
  U_update_copy=copy.copy(U_update)
  U_first=U_update
  E2=0
  for i in xrange(Max_SVD_iteratoin):
   E1=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN)
   #print E1, i, '\n'
   if E1>E2 or i is 0:
    U_update_copy=copy.copy(U_update)
    E2=E1
   else:
    U_update=U_update_copy
    print 'Notoptimized=', i
    U_list[L_position][L_lay_selected].putBlock(U_update)
    break
   svd=Env_Uni_inner[L_position][L_lay_selected].transpose().getBlock().svd()
   temporary_matrix=svd[0]*svd[2]
   U_update=temporary_matrix.transpose()
   U_list[L_position][L_lay_selected].putBlock(U_update)
  Energy_f=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN)
  print 'E_s=', Energy_s, '\n', 'E_f=', Energy_f
  if Energy_s > Energy_f:
   print 'Notoptimized= E > E1',  Energy_s, Energy_f
   U_list[L_position][L_lay_selected].putBlock(U_first)

 elif Method is 'SteepestDescent':
  E=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN)
  #print 'E=', E
  
  E2=0.0
  U_update=copy.copy(U_list[L_position][L_lay_selected].getBlock())
  U_first=copy.copy(U_update)
  Gamma= 1.00
  count=0
  for i in xrange(Max_Steepest_iteratoin):
   U_list[L_position][L_lay_selected].putBlock(U_update)
   E1=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN)
   #print 'E=', E1, count
   count=0
   #print 'test', (Env_Uni_inner[L_position][L_lay_selected].transpose().getBlock()*U_update).trace()
   D_u=-4.00*Env_Uni_inner[L_position][L_lay_selected].getBlock()
   D_u_trans=copy.copy(D_u)
   D_u_trans.transpose()
   #print D_u, D_u_trans, (-1.00)*D_u, U_update*D_u_trans*U_update 
   Z_decent=U_update*D_u_trans*U_update+(-1.00)*D_u
   Z_decent_trans=copy.copy(Z_decent)
   Z_decent_trans.transpose()
   Norm_Z=(Z_decent_trans*Z_decent).trace() / 2.00
   if Norm_Z < 1.0e-7:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   while Break_loop is 1:
    count+=1
    Temporary=U_update+(2.00)*Gamma*Z_decent
    svd=Temporary.svd()
    Temporary=svd[0]*svd[2]
    U_list[L_position][L_lay_selected].putBlock(Temporary)
    #print U_list[L_position][L_lay_selected], Temporary
    E2=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN)
    #print 'E2=', E2, E1-E2, -Norm_Z*Gamma 
    if E1-E2 <=(-Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0
   
   Break_loop=1
   while Break_loop is 1:
    count+=1
    #print 'Numbers=', count 
    Temporary=U_update+Gamma*Z_decent
    svd=Temporary.svd()
    Temporary=svd[0]*svd[2]
    U_list[L_position][L_lay_selected].putBlock(Temporary)
    E2=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN)
    #print 'E2s=', E2, 'Gamma=', Gamma
    #print 'abs=', abs((0.5)*Norm_Z*Gamma) 
    if abs((0.5)*Norm_Z*Gamma) <1.0e-11 or  abs(E1-E2)<1.0e-11 :
     #print 'break, Gamma is too small', 'E1-E2=', abs(E1-E2)
     Temporary=U_update+Gamma*Z_decent
     svd=Temporary.svd()
     Temporary=svd[0]*svd[2]
     break
    #print E1-E2, (-0.5)*Norm_Z*Gamma  
    if E1-E2 > (-0.5)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0

   Temporary=U_update+Gamma*Z_decent
   svd=Temporary.svd()
   U_update=svd[0]*svd[2]
  
  U_list[L_position][L_lay_selected].putBlock(U_update)
  E_f=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN)
  
  #print 'E_s=', Energy_s, '\n', 'E_f=', E_f 
  if Energy_s > E_f:
   print 'Notoptimized= E > E1', E_f,  Energy_s
   U_list[L_position][L_lay_selected].putBlock(U_first)


    
    
    
  
  
 
 
