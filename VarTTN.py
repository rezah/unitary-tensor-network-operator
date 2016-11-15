import pyUni10 as uni10
#import syss
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time
import mpo
import env
import optimize 
import plotdata
################################## Set parameter #######################################
Tech='Regular'                  # MERA, Regular  
Model='Heisenberg'              #could be: 'Ising' and or 'Heisenberg'
L=16                            #Numbers of particles should be even!
L_lay=[0,1,2,3,4]               #Numbers of layers <= 5
realizations=50                  #Number of realizations        
Accuracy=1.00e-7                #Accuracy of variational variance  
d=2                              #pysical bond-dimension  
chi=5                            #bond-dimension of MPO
W=8                              #random interval, [-W,W]
J=1.0                            #coupling, 1.0
J2=0.30                          #coupling, 0.3
Fieldz=0.60                      #Field in z direction, 0.6
hz_list=[]                       # list of randomness 
U_delta=0.020                     #if it's zero, U_list is intialized by Identity
Method='SteepestDescent'         #methods: CGarmjo, CGpoly, SVD, SteepestDescent, SteepestDescentploy 
Randomness='Not-Fixed'               #Fixed
Max_number_iteratoin_SVD=400       # maximum number of sweeps for SVD method 
Max_number_iteratoin_Steepest=400  # maximum number of sweeps for SteepestDescent method
Max_number_iteratoin_CG=400         #maximum number of sweeps for CG method

Max_SVD_iteratoin=16               #maximum number of SVD iteration
Max_Steepest_iteratoin=16          #maximum number of SteepestDescent iteration
Max_CG_iteratoin=16                #maximum number of SteepestDescent iteration
#######################################################################################3
def Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list, Gamma,Tech,trH2):

 count=0
 Loop_iter=0
 E0=0

 while Loop_iter is 0: 
  for j in xrange(1):
   count+=1
   for i in xrange( L/2 ):
    #make Env right
    if i is 0:
     for m in xrange((L/2)-1):
      L_position1=(L/2)-1-m
      #print L_position
      Environment_Right[L_position1]=env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position1, d, Environment_Right)

    #print i
    L_position=i
    optimize.optimize_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list,Gamma,Tech)
    if i is not L/2 -1:
     Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
  E1=E_list[len(E_list)-1]
  if abs((E0-E1)/E1) < Accuracy:Loop_iter=1;
  if count > Max_number_iteratoin_SVD : Loop_iter=1;
  print 'count=', count, 'Accuracy=', abs((E0-E1)/E1),'E0=', E0,'E1=', E1, 'Loop_iter=', Loop_iter, 'Count_list=', Count_list[len(Count_list)-1], 'Variance=', (trH2-E1) * (1.0/(2**L))
  E0=E1
########################################################################################
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

#print hz_list

##########    only works for Heisenberg    ##########################
#trH2=mpo.print_trH2(L, hz_list,J)
#print trH2 
#mpo.avarage_Energy_power_2(L, hz_list,J)

######### retrun MPO list, 1: one-site MPO, 2: two-site MPO; mpo_uni10_bl, mpo_uni10_br, mpo_uni10_bl2, mpo_uni10_br2###########################################
mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model )
trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L )
print trH2

U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta)
copy_U_list=optimize.copy_U_list_function(U_list, L, L_lay)

######### make a mpo representation of unitary U, up stand for U, down for U^{T} ######### 
mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down', Tech)


perl_label_up=[None]*4
perl_label_up=[None]*4
Bond_IN=[None]*4
List_position=[ 0, 1, 2, (L/2)-1]
for i in xrange(4):
 L_position=List_position[i]
 perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)




Environment_Left=[None]*(L/2)
Environment_Right=[None]*(L/2)
Environment_Uni=[None]*(L/2)
Env_Uni_inner, Gamma=Initialize_function(L, L_lay)



file = open("Data/varianceAll.txt", "w")
variance_final1=[]
variance_final2=[]
Count_final2=[]
Count_final1=[]
for q in xrange(realizations):
 print '\n', 'q_regular=', q, '\n'
 variance_list1=[]
 variance_list2=[]

 E_list1=[]
 Count_list1=[]
 E_list2=[]
 Count_list2=[]
 Randomness='Not-Fixed'
 Tech='Regular'
 L_lay=[0,1,2,3,4]
 mpo.intialize_coupling(L, hz_list, W, Randomness, Model)
 mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model )
 trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L )

 Method='SVD'
 U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta)
 mpo_U_list_up = mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
 mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down', Tech)
 Env_Uni_inner, Gamma=Initialize_function(L, L_lay)
 for i in xrange(4):
  L_position=List_position[i]
  perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)
 E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
 Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list1, Count_list1,Gamma,Tech,trH2 )
 

 ###################################   MERA  ##########################
 print '\n', 'q_MERA=', q, '\n'
 Tech='MERA'
 L_lay=[0,1,2,3,4,5,6,7]
 for i in xrange(4):
  L_position=List_position[i]
  perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)
 Env_Uni_inner, Gamma=Initialize_function(L, L_lay)
 U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta)
 copy_U_list1=optimize.copy_U_list_function(U_list, L, L_lay)
 U_list=copy_U_list1

 mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
 mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down', Tech)
 E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)


 Method='SVD'

 Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list2, Count_list2,Gamma,Tech,trH2 )

 
 variance_list1=[ (trH2-E_list1[i])*(1.0/(2**L)) for i in xrange(len(E_list1)) ]
 variance_list2=[ (trH2-E_list2[i])*(1.0/(2**L)) for i in xrange(len(E_list2)) ]
 variance_final1.append(variance_list1[len(variance_list1)-1])
 variance_final2.append(variance_list2[len(variance_list2)-1])
 Count_final1.append(Count_list1[len(Count_list1)-1])
 Count_final2.append(Count_list2[len(Count_list2)-1])

 plotdata.plot(Count_list1, variance_list1,Count_list2, variance_list2, q,L,hz_list)
 plotdata.Store(variance_final1, variance_final2, q , L , file, Count_final1, Count_final2)
 
#############################################################################################################33

file.close()


