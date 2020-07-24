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
import numpy as np
from scipy.linalg import expm
from numpy.linalg import inv
N_poly=3
q_order=4


def Mat_np_to_Uni(Mat_np):
 d0=np.size(Mat_np,0)
 d1=np.size(Mat_np,1)
 Mat_uni=uni10.Matrix(d0,d1)
 for i in xrange(d0):
  for j in xrange(d1):
   Mat_uni[i*d1+j]=Mat_np[i,j]
 return  Mat_uni
 
def Mat_uni_to_np(Mat_uni):
 dim0=int(Mat_uni.row())
 dim1=int(Mat_uni.col())
 Mat_np=np.zeros((dim0,dim1))
 for i in xrange(dim0):
  for j in xrange(dim1):
   Mat_np[i,j]=Mat_uni[i*dim1+j]
 return  Mat_np

##################  Line-Poly ##############################################
def Line_search_pol(H_direct, U_update,count,U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech):
 ##print H_direct
 

 dim0=int(H_direct.row())
 dim1=int(H_direct.col())
 Hk=np.zeros((dim0,dim1))
 for i in xrange(dim0):
  for j in xrange(dim1):
   Hk[i,j]=H_direct[i*dim1+j]
#  
 dim0=int(U_update.row())
 dim1=int(U_update.col())
 Wk=np.zeros((dim0,dim1))
 for i in xrange(dim0):
  for j in xrange(dim1):
   Wk[i,j]=U_update[i*dim1+j]
   


 A=np.linalg.eig(Hk)
 eigs=abs(A[0])
 maxeig=np.amax(eigs)
 ##print 'maxeig', maxeig, eigs
 T_mu=2.00*np.pi/(q_order*maxeig);
 ##print 'q_order=',q_order, N_poly
 ##print 'T_mu=', T_mu 
 
 mu_step_poly=T_mu/N_poly
 R_poly=expm(-mu_step_poly*Hk)
 R_mu_poly=np.eye(np.size(Wk,0))
 d1_poly=np.zeros((N_poly+1,1))
 
# #print 'np.size(Wk,0)=', np.size(Wk,0)
# Temporary=Mat_np_to_Uni(np.dot(R_mu_poly,Wk))
# #print 'Wk', Wk,'\n',np.dot(R_mu_poly,Wk), '\n', Temporary
# U_list[L_position][L_lay_selected].putBlock(Temporary)
# #print U_list[L_position][L_lay_selected]
# E2=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
# #print 'E2=', E2, (U_list[L_position][L_lay_selected]*Env_Uni_inner[L_position][L_lay_selected])[0]
# Temporary=4.00*Env_Uni_inner[L_position][L_lay_selected].getBlock(); Dk=Mat_uni_to_np(Temporary)
# #print Temporary, Dk
# 
 for n_poly  in  xrange(0,N_poly+1):
  Mat_np1=np.dot(R_mu_poly,Wk) 
  Temporary=Mat_np_to_Uni(Mat_np1)
  U_list[L_position][L_lay_selected].putBlock(Temporary)
  E2=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
  Temporary=4.00*Env_Uni_inner[L_position][L_lay_selected].getBlock()
  Dk=Mat_uni_to_np(Temporary)
  d1_poly[n_poly,0]=2.00*np.real(np.trace(np.dot(Dk,np.dot(np.dot(np.transpose(Wk),np.transpose(R_mu_poly)),np.transpose(Hk)))))
  R_mu_poly=np.dot(R_mu_poly, R_poly)


 ##print 'd1_poly', d1_poly
 mu_poly=np.transpose(np.linspace(0, T_mu, N_poly+1)) 
 C=np.zeros((N_poly,N_poly))
 
 for i  in xrange(N_poly):
  for j  in xrange(N_poly):
   C[j,i]= mu_poly[j+1]**(i+1)
 ##print 'C=', C,'\n','\n','\n','\n'

 ##print inv(C),'\n','\n','\n','\n'


 Array=np.zeros((N_poly, 1))
 for i in xrange(N_poly):
  Array[i]=d1_poly[i+1]-d1_poly[0]
 
 a=np.dot(inv(C) , Array)
 a_1=np.append(d1_poly[0], a)
 a_2=np.reshape(a_1,[N_poly+1])
 
 ##print 'a=', a_2,'\n','\n','\n'
 ##print Array
 ##print np.flipud(a_2)
 Results=np.roots(np.flipud(a_2))
 ##print Results
 Real_roots=Results[np.isreal(Results)]
 Real_roots=np.real(Real_roots)
 Real_positive=Real_roots[Real_roots>0]
 Real_positive=np.amin(Real_positive)
 
 if Real_positive < 1.0e-12:  
  Real_positive=0; #print 'Root is very small';
 count=count+(N_poly+1)
 return Real_positive , count


#################################   Line Search    #################################
def Line_search(Z_decent, Gamma, E1, U_update,count, Norm_Z,U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech):
 Break_loop=1
 Temporary=exp_matrix(Z_decent,-Gamma,U_update)
 while Break_loop is 1:
  count+=1
  Temporary=Temporary*Temporary
  Temporary_trans=copy.copy(Temporary)
  U_list[L_position][L_lay_selected].putBlock(Temporary)
  E2=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
  if E1-E2 <= -Norm_Z*Gamma:
   Gamma*=1.000
  else:
   Break_loop=0
 Break_loop=1

 while Break_loop is 1:
  count+=1
  Temporary=exp_matrix(Z_decent,-Gamma, U_update)
  Temporary_trans=copy.copy(Temporary)
  U_list[L_position][L_lay_selected].putBlock(Temporary)
  E2=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
  if Gamma < 1.0e-11:
   Break_loop=0
   
  if E1-E2 > (-0.5)*Norm_Z*Gamma:
   Gamma*=0.50
  else:
   Break_loop=0
 return Gamma , count
############################################Exponential################################################
def exp_matrix(H_direct,Gamma, U_update):
 dim0=int(H_direct.row())
 dim1=int(H_direct.col())
 M0=np.zeros((dim0,dim1))
 M_return1=uni10.Matrix(dim0,dim1)
 #t0 = time.time()
 for i in xrange(dim0):
  for j in xrange(dim1):
   M0[i,j]=H_direct[i*dim1+j]
 ##print time.time() - t0, "replacing"

 M0=expm(+Gamma*M0)
 
 for i in xrange(dim0):
  for j in xrange(dim1):
   M_return1[i*dim1+j]=M0[i,j]
 M_return1=M_return1*U_update
 return M_return1 
################################      Energy           ##########################
def Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech):
 
 Update_Unitary_change(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L,Tech)
 
 Environment_Uni[L_position]=env.Environment_uni_function(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right)
 
 env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner,Tech )
 E=(Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0]
 return E

####################     Update mpo_U_up/down & Environment_Left    ##################
def Update_Unitary_change(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L,Tech):
 mpo_U_up=mpo.make_mpo_U(U_list, L_position, L_lay, L, 'up', Tech)
 mpo_U_down=mpo.make_mpo_U(U_list, L_position, L_lay, L, 'down', Tech)
 mpo_U_list_up[L_position]=mpo_U_up
 mpo_U_list_down[L_position]=mpo_U_down
# Environment_Left[L_position]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
################################   Copy U_list       ###################################
def copy_U_list_function(U_list,L,L_lay,Tech):

 if (Tech is 'MERA') and (len(L_lay) is 6):
  L_m=L/2
 elif (Tech is 'MERA') and (len(L_lay) is 9):
  L_m=((L-8)/3)+4
 else: L_m=L;
 
 
 U_list_copy=[]
 
 for i in xrange(L_m/2):
  U_list_copy.append([])
 for i in xrange(L_m/2):
  for j in xrange(len(L_lay)):
    U_list_copy[i].append([None])
 for i in xrange(L_m/2):
   for j in xrange(len(L_lay)):
    U_list_copy[i][j]=copy.copy(U_list[i][j])
 return U_list_copy
#############################################################################################

###################################Optimze  ###############################################
def optimize_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list,Gamma,Tech):

 if Tech is 'Regular':
  if L_position is not L/2-1:
   for i in xrange(len(L_lay)):
    L_lay_selected=i
    ##print 'L_position, L, L_lay', L_position, L, L_lay_selected
    optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma,Tech)
  elif L_position is L/2 -1:
   for i in xrange(0,len(L_lay),2):
    L_lay_selected=i
    ##print 'L_position, L, L_lay', L_position, L, L_lay_selected
    optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma,Tech)
 elif Tech is 'MERA':
  if (L_position is 0) and (len(L_lay) is 8) :
   for i in xrange(len(L_lay)-3):
    L_lay_selected=i
    ##print 'L_position, L, L_lay', L_position, L, L_lay_selected
    optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma,Tech)
  elif  (L_position is 0) and ( (len(L_lay) is 6) or (len(L_lay) is 9) ) :
   for i in xrange(2):
    L_lay_selected=i
    ##print 'L_position, L, L_lay', L_position, L, L_lay_selected
    optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma,Tech)
  
  

  if (L_position is not  0) and (L_position is not L/2 -1) and (len(L_lay) is 8) :
   if L_position % 2 is 0:
    for i in xrange(len(L_lay)):
     L_lay_selected=i
     ##print 'L_position, L, L_lay', L_position, L, L_lay_selected
     optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma,Tech)
   if L_position % 2 is 1:
    for i in xrange(len(L_lay)-5):
     L_lay_selected=i
     ##print 'L_position, L, L_lay', L_position, L, L_lay_selected
     optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma,Tech)
  elif (L_position is not  0) and (L_position is not L/2 -1) and (len(L_lay) is 6) :
   for i in xrange(len(L_lay)):
     L_lay_selected=i
     ##print 'L_position, L, L_lay', L_position, L, L_lay_selected
     optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma,Tech)
  elif (L_position is not  0) and (L_position is not L/2 -1) and (len(L_lay) is 9) :
   for i in xrange(len(L_lay)-1):
     L_lay_selected=i
     ##print 'L_position, L, L_lay', L_position, L, L_lay_selected
     optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma,Tech)
  
  if L_position is L/2 -1 and (len(L_lay) is 8) :
   for i in xrange(0,(len(L_lay)-3),2):
    L_lay_selected=i
    ##print 'L_position, L, L_lay', L_position, L, L_lay_selected
    optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma,Tech)
  elif (L_position is L/2 -1) and ( (len(L_lay) is 6) or (len(L_lay) is 9) )  :
   for i in xrange(3):
    L_lay_selected=i
    ##print 'L_position, L, L_lay', L_position, L, L_lay_selected
    optimize_inner_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma,Tech)


################################ Optimize inner function #########################
def optimize_inner_function(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position, L_lay_selected,Method,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin,E_list, Count_list,Gamma_list,Tech):
 Energy_s=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
 ###print 'E_s=', Energy_s
 
 if Method is 'SVD':
  U_update=copy.copy(U_list[L_position][L_lay_selected]).getBlock()
  U_update_copy=copy.copy(U_update)
  U_first=U_update
  E2=0
  
  if len(Count_list) is 0: count=0; 
  else: count = Count_list[len(Count_list)-1];
  
  for i in xrange(Max_SVD_iteratoin):
   count+=1
   E1=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
   if E1>E2 or i is 0:
    U_update_copy=copy.copy(U_update)
    if abs((E2-E1)/E1) < 1.0e-9:
     ###print E2, E1, abs((E2-E1)/E1), i
     break
    E2=E1
   else:
    U_update=U_update_copy
    ##print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
    U_list[L_position][L_lay_selected].putBlock(U_update)
    break
   svd=Env_Uni_inner[L_position][L_lay_selected].transpose().getBlock().svd()
   temporary_matrix=svd[0]*svd[2]
   U_update=temporary_matrix.transpose()
   U_list[L_position][L_lay_selected].putBlock(U_update)
  Energy_f=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
  Count_list.append(count)
  E_list.append(Energy_f)
  #print 'E_s=', Energy_s, '\n', 'E_f=', Energy_f
  if Energy_s > Energy_f:
   print 'Notoptimized= E > E1, SVD',  Energy_s, Energy_f
   U_list[L_position][L_lay_selected].putBlock(U_first)

 elif (Method is 'CGarmjo') or (Method is  'CGpoly') or (Method is  'SteepestDescentploy'):
  E=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
  ###print 'E=', E
  E2=0.0
  U_update=copy.copy(U_list[L_position][L_lay_selected].getBlock())
  U_first=copy.copy(U_update)
  Gamma= Gamma_list[L_position][L_lay_selected]
  #Gamma=1.00
  if len(Count_list) is 0: count=0;
  else: count = Count_list[len(Count_list)-1];
  E_previous=0
  ###print 'Max_CG_iteratoin', Max_CG_iteratoin 
  for i in xrange(Max_CG_iteratoin):
   if (i % (d*d*d*d)) is 0:
    count+=1
    U_list[L_position][L_lay_selected].putBlock(U_update)
    E1=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
    #E_list.append(E1)
    #Count_list.append(count)
    D_u=4.00*Env_Uni_inner[L_position][L_lay_selected].getBlock()
    D_u_trans=copy.copy(D_u)
    D_u_trans.transpose()
    U_update_trans=copy.copy(U_update)
    U_update_trans.transpose() 
    Z_decent=(-1.00)*(D_u*U_update_trans+(-1.00)*U_update*D_u_trans)
    H_direct=copy.copy(Z_decent)
   Z_decent_trans=copy.copy(Z_decent)
   Z_decent_trans.transpose()
   H_direct_trans=copy.copy(H_direct)
   H_direct_trans.transpose()
   Norm_Z=Z_decent_trans*Z_decent
   Norm_Z=Norm_Z.trace() / 2.00
   if Norm_Z < 1.0e-7:
    ###print 'Break Norm=', Norm_Z
    break
   Norm_Z=H_direct_trans*Z_decent
   Norm_Z=Norm_Z.trace() / 2.00
   Gamma=1.00
   if Method is 'CGarmjo':
    Gamma , count=Line_search(H_direct, Gamma, E1, U_update,count, Norm_Z,U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
   elif (Method is 'CGpoly') or (Method is  'SteepestDescentploy'):
    Gamma , count=Line_search_pol(H_direct, U_update,count,U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)


   ###print 'Gamma=', Gamma, count
   if E1>E_previous or i is 0:
    if abs((E_previous-E1)/E1) < 1.0e-8:
     ###print E_previous, E1, abs((E_previous-E1)/E1), i
     break
    E_previous=E1



   Temporary=exp_matrix(H_direct,-Gamma, U_update)
   U_update=Temporary
   U_update_trans=copy.copy(U_update)
   U_update_trans.transpose()
   U_list[L_position][L_lay_selected].putBlock(U_update)
   E1=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
   D_u=4.00*Env_Uni_inner[L_position][L_lay_selected].getBlock()
   D_u_trans=copy.copy(D_u)
   D_u_trans.transpose()
   Z_decent1=(-1.00)*(D_u*U_update_trans+(-1.00)*U_update*D_u_trans)
   Z_decent1_trans=copy.copy(Z_decent1)
   Z_decent1_trans.transpose()
   Z_decent_h=Z_decent1+(-1.00)*Z_decent
   A=(Z_decent_h*Z_decent1_trans).trace()
   B=(Z_decent*Z_decent_trans).trace()
   norm_gamma=A/B
   #norm_gamma=0
   ###print norm_gamma 
   if Method is  'SteepestDescentploy': norm_gamma=0.0;
   H_direct1=Z_decent1+(1.0)*norm_gamma*H_direct
   Check=(H_direct1*Z_decent1_trans).trace()
   if Check<0: H_direct1=Z_decent1;
   H_direct=copy.copy(H_direct1)
   Z_decent=copy.copy(Z_decent1)
   count+=1
   #E_list.append(E1)
   #Count_list.append(count)

  U_list[L_position][L_lay_selected].putBlock(U_update)
  E_f=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
  Gamma_list[L_position][L_lay_selected]=Gamma+0.01*Gamma
  E_list.append(E_f)
  Count_list.append(count)

  #print 'E_s=', Energy_s, '\n', 'E_f=', E_f 
  if Energy_s > E_f:
   #print 'Notoptimized= E > E1, CGPOLY', E_f,  Energy_s
   U_list[L_position][L_lay_selected].putBlock(U_first)
 
 
 elif Method is 'SteepestDescent':
  E=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
  ##print 'E=', E
  
  E2=0.0
  U_update=copy.copy(U_list[L_position][L_lay_selected].getBlock())
  U_first=copy.copy(U_update)
  Gamma= Gamma_list[L_position][L_lay_selected]
  #Gamma=1.0
  if len(Count_list) is 0: count=0 
  else: count = Count_list[len(Count_list)-1]
  E_previous=0
  for i in xrange(Max_Steepest_iteratoin):
   count+=1
   U_list[L_position][L_lay_selected].putBlock(U_update)
   E1=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
   ##print 'E=', E1, count
   #E_list.append(E1)
   #Count_list.append(count)
   ##print 'test', (Env_Uni_inner[L_position][L_lay_selected].transpose().getBlock()*U_update).trace()
   D_u=-4.00*Env_Uni_inner[L_position][L_lay_selected].getBlock()
   D_u_trans=copy.copy(D_u)
   D_u_trans.transpose()
   ##print D_u, D_u_trans, (-1.00)*D_u, U_update*D_u_trans*U_update 
   Z_decent=U_update*D_u_trans*U_update+(-1.00)*D_u
   Z_decent_trans=copy.copy(Z_decent)
   Z_decent_trans.transpose()
   Norm_Z=(Z_decent_trans*Z_decent).trace() / 2.00
    
    
   if E1>E_previous or i is 0:
    if abs((E_previous-E1)/E1) < 1.0e-8:
     ##print E_previous, E1, abs((E_previous-E1)/E1), i
     break
   E_previous=E1
   
   if Norm_Z < 1.0e-7:
    ##print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   Gamma=1.0
   while Break_loop is 1:
    count+=1
    Temporary=U_update+(2.00)*Gamma*Z_decent
    svd=Temporary.svd()
    Temporary=svd[0]*svd[2]
    U_list[L_position][L_lay_selected].putBlock(Temporary)
    ##print U_list[L_position][L_lay_selected], Temporary
    E2=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
    ##print 'E2=', E2, E1-E2, -Norm_Z*Gamma 
    if E1-E2 <=(-Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0
   
   Break_loop=1
   while Break_loop is 1:
    count+=1
    ##print 'Numbers=', count 
    Temporary=U_update+Gamma*Z_decent
    svd=Temporary.svd()
    Temporary=svd[0]*svd[2]
    U_list[L_position][L_lay_selected].putBlock(Temporary)
    E2=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
    ##print 'E2s=', E2, 'Gamma=', Gamma
    ##print 'abs=', abs((0.5)*Norm_Z*Gamma) 
    if abs((0.5)*Norm_Z*Gamma) <1.0e-11 or  abs(E1-E2)<1.0e-11 :
     ##print 'break, Gamma is too small', 'E1-E2=', abs(E1-E2)
     Temporary=U_update+Gamma*Z_decent
     svd=Temporary.svd()
     Temporary=svd[0]*svd[2]
     break
    ##print E1-E2, (-0.5)*Norm_Z*Gamma  
    if E1-E2 > (-0.5)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0

   Temporary=U_update+Gamma*Z_decent
   svd=Temporary.svd()
   U_update=svd[0]*svd[2]
  
  U_list[L_position][L_lay_selected].putBlock(U_update)
  E_f=Energy_newunitary(U_list, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, L_lay,L, Environment_Right,L_lay_selected,Env_Uni_inner,Environment_Uni, perl_label_up, Bond_IN,Tech)
  Gamma_list[L_position][L_lay_selected]=Gamma+0.01*Gamma
  E_list.append(E_f)
  Count_list.append(count)
  ##print 'E_s=', Energy_s, '\n', 'E_f=', E_f 
  if Energy_s > E_f:
   #print 'Notoptimized= E > E1, Steepest', E_f,  Energy_s
   U_list[L_position][L_lay_selected].putBlock(U_first)


    
    
    
  
  
 
 
