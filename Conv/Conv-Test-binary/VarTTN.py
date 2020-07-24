import pyUni10 as uni10
#import sys
#import numpy as np
#import matplotlib
#matplotlib.use('pdf')
#import matplotlib.pyplot as plt
#import pylab
import random
import copy
import time
import mpo
import env
import optimize 
import plotdata
################################## Set parameter #######################################
Tech='MERA'                  # MERA, Regular  
Model='Heisenberg'              #could be: 'Ising' and or 'Heisenberg'
L=32                     #Numbers of particles should be 4*k and also k={1,2,3,4,...}!6*(2k)+8
L_lay=[0,1,2,3,4,5]               #Numbers of layers <= 5
realizations=30                  #Number of realizations        
Accuracy=8.00e-6                #Accuracy of variational variance  
d=2                              #pysical bond-dimension  
chi=5                            #bond-dimension of MPO
W=6                              #random interval, [-W,W]
J=1.0                            #coupling, 1.0
J2=0.30                          #coupling, 0.3
Fieldz=0.60                      #Field in z direction, 0.6
hz_list=[]                       # list of randomness 
U_delta=0.200                     #if it's zero, U_list is intialized by Identity
Method='SVD'         #methods: CGarmjo, CGpoly, SVD, SteepestDescent, SteepestDescentploy 
Randomness='Non-Fixed'               #Fixed
Max_number_iteratoin_SVD=22       # maximum number of sweeps for SVD method 
Max_number_iteratoin_Steepest=22  # maximum number of sweeps for SteepestDescent method
Max_number_iteratoin_CG=22         #maximum number of sweeps for CG method

Max_SVD_iteratoin=35               #maximum number of SVD iteration
Max_Steepest_iteratoin=25          #maximum number of SteepestDescent iteration
Max_CG_iteratoin=25                #maximum number of SteepestDescent iteration
#######################################################################################3
L_fixed=L
def Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list, Gamma,Tech,trH2,E_time1, time1 ):

 count=0
 Loop_iter=0
 E0=0
 t0=time.time()
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
    E_time1.append(   (trH2-E_list[len(E_list)-1])*((1.0/(2**(L_fixed))))    )
    time1.append(  (time.time()- t0) / 10000  )

    if i is not L/2 -1:
     Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
  
  E1=E_list[len(E_list)-1]
  if abs((E0-E1)/E1) < Accuracy:Loop_iter=1;
  if count > Max_number_iteratoin_SVD : Loop_iter=1;
  print 'count=', count, 'Accuracy=', abs((E0-E1)/E1),'E0=', E0,'E1=', E1, trH2, 'Loop_iter=', Loop_iter, 'Count_list=', Count_list[len(Count_list)-1], 'Variance=', (trH2-E1) * (1.0/(2**(L_fixed)))
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
trH2=mpo.print_trH2(L, hz_list,J,Model)
#print 'trH2:=', trH2 
variance_Energy=mpo.avarage_Energy_power_2(L, hz_list,J)
#print variance_Energy
######### retrun MPO list, 1: one-site MPO, 2: two-site MPO; mpo_uni10_bl, mpo_uni10_br, mpo_uni10_bl2, mpo_uni10_br2###########################################
mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2, Model, Tech, L_lay )

#print mpo_list2[0].printDiagram(),mpo_list2[1].printDiagram(),mpo_list2[2].printDiagram(),mpo_list2[3].printDiagram(),mpo_list2[4].printDiagram(),mpo_list2[5].printDiagram()

trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L, Tech, L_lay )
#print 'trH2:=', trH2

U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta,Tech)

#print U_list[0][1].printDiagram(), U_list[0][0].printDiagram(), U_list[((L-8)/6)+1][0].printDiagram(), U_list[((L-8)/6)+1][1].printDiagram(),U_list[((L-8)/6)+1][2].printDiagram(), U_list[1][1].printDiagram(), U_list[1][6].printDiagram(), U_list[1][7].printDiagram()


copy_U_list=optimize.copy_U_list_function(U_list, L, L_lay,Tech)

######### make a mpo representation of unitary U, up stand for U, down for U^{T} ######### 



perl_label_up=[None]*4
perl_label_up=[None]*4
Bond_IN=[None]*4


if (Tech is 'MERA') and (len(L_lay) is 6):
 L_m=L/2
elif (Tech is 'MERA') and (len(L_lay) is 9):
 L_m=((L-8)/3)+4
else: L_m=L;


mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L_m, 'up', Tech)
mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L_m, 'down', Tech)

List_position=[ 0, 1, 2, (L_m/2)-1]
for i in xrange(4):
 L_position=List_position[i]
 perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)



Environment_Left=[None]*(L_m/2)
Environment_Right=[None]*(L_m/2)
Environment_Uni=[None]*(L_m/2)
Env_Uni_inner, Gamma=Initialize_function(L_m, L_lay)


file = open("Data/varianceAll.txt", "w")
variance_final1=[]
variance_final2=[]
variance_final3=[]
variance_final4=[]

Count_final1=[]
Count_final2=[]
Count_final3=[]
Count_final4=[]

time_final1=[]
time_final2=[]
time_final3=[]
time_final4=[]



for q in xrange(realizations):
 print '\n', 'MERA_binary=', q, '\n'
 variance_list1=[]
 variance_list2=[]

 E_list1=[]
 Count_list1=[]
 E_list2=[]
 Count_list2=[]
 E_list3=[]
 Count_list3=[]
 E_list4=[]
 Count_list4=[]

 
 E_time1=[]
 time1=[]
 E_time2=[]
 time2=[]
 E_time3=[]
 time3=[]
 E_time4=[]
 time4=[]
 
 
 #Randomness='Non-Fixed'
 mpo.intialize_coupling(L, hz_list, W, Randomness, Model)
 print hz_list
 trH2=mpo.print_trH2(L, hz_list,J, Model)
 print 'trH2=', trH2 
 
 
 
 Tech='MERA'
 #L_lay=[0,1,2,3,4,5]
 mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model,Tech,L_lay )
 trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L, Tech,L_lay)
 print '\n', 'trH2=', trH2
 Method='SVD'
 #Method='CGpoly'
 print Method
 U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta,Tech)
 U_list_copy=optimize.copy_U_list_function(U_list,L,L_lay,Tech)
 mpo_U_list_up = mpo.make_mpo_U_list(U_list, L_lay, L/2, 'up', Tech)
 mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L/2, 'down', Tech)

 #print mpo_U_list_up[0].printDiagram(), U_list[0][0].printDiagram(),U_list[0][1].printDiagram()

 List_position=[ 0, 1, 2, (L/4)-1]
 for i in xrange(4):
  L_position=List_position[i]
  perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)

 Environment_Left=[None]*(L/4)
 Environment_Right=[None]*(L/4)
 Environment_Uni=[None]*(L/4)
 Env_Uni_inner, Gamma=Initialize_function(L/2, L_lay)
 E_f=mpo.Energy_through_env(U_list, L_lay, L/2, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
 print 'E_0 = ', E_f,'\n'
 
 Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L/2,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list1, Count_list1,Gamma,Tech,trH2,E_time1, time1)
 
###################

 print hz_list
 trH2=mpo.print_trH2(L, hz_list,J, Model)
 print 'trH2=', trH2 
 
 
 
 Tech='MERA'
 #L_lay=[0,1,2,3,4,5]
 mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model,Tech,L_lay )
 trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L, Tech,L_lay)
 print '\n', 'trH2=', trH2
# Method='SVD'
 Method='CGpoly'
 print Method
 #U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta,Tech)
 U_list=optimize.copy_U_list_function(U_list_copy,L,L_lay,Tech)

 mpo_U_list_up = mpo.make_mpo_U_list(U_list, L_lay, L/2, 'up', Tech)
 mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L/2, 'down', Tech)

 #print mpo_U_list_up[0].printDiagram(), U_list[0][0].printDiagram(),U_list[0][1].printDiagram()

 List_position=[ 0, 1, 2, (L/4)-1]
 for i in xrange(4):
  L_position=List_position[i]
  perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)

 Environment_Left=[None]*(L/4)
 Environment_Right=[None]*(L/4)
 Environment_Uni=[None]*(L/4)
 Env_Uni_inner, Gamma=Initialize_function(L/2, L_lay)
 E_f=mpo.Energy_through_env(U_list, L_lay, L/2, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
 print 'E_0 = ', E_f,'\n'
 
 Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L/2,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list2, Count_list2,Gamma,Tech,trH2,E_time2, time2)
##################




 print hz_list
 trH2=mpo.print_trH2(L, hz_list,J, Model)
 print 'trH2=', trH2 
 
 
 
 Tech='MERA'
 #L_lay=[0,1,2,3,4,5]
 mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model,Tech,L_lay )
 trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L, Tech,L_lay)
 print '\n', 'trH2=', trH2
# Method='SVD'
 Method='CGarmjo'
 print Method
 #U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta,Tech)
 U_list=optimize.copy_U_list_function(U_list_copy,L,L_lay,Tech)

 mpo_U_list_up = mpo.make_mpo_U_list(U_list, L_lay, L/2, 'up', Tech)
 mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L/2, 'down', Tech)

 #print mpo_U_list_up[0].printDiagram(), U_list[0][0].printDiagram(),U_list[0][1].printDiagram()

 List_position=[ 0, 1, 2, (L/4)-1]
 for i in xrange(4):
  L_position=List_position[i]
  perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)

 Environment_Left=[None]*(L/4)
 Environment_Right=[None]*(L/4)
 Environment_Uni=[None]*(L/4)
 Env_Uni_inner, Gamma=Initialize_function(L/2, L_lay)
 E_f=mpo.Energy_through_env(U_list, L_lay, L/2, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
 print 'E_0 = ', E_f,'\n'
 
 Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L/2,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list3, Count_list3,Gamma,Tech,trH2,E_time3, time3)
##################

 
 variance_list1=[ (trH2-E_list1[i])*(1.0/(2**L)) for i in xrange(len(E_list1)) ]
 variance_list2=[ (trH2-E_list2[i])*(1.0/(2**L)) for i in xrange(len(E_list2)) ]
 variance_list3=[ (trH2-E_list3[i])*(1.0/(2**L)) for i in xrange(len(E_list3)) ]
#  variance_list4=[ 1 ]
# variance_list2=[ 1 ]
# variance_list3=[ 1 ]
 variance_list4=[ 1 ]

 variance_final1.append(variance_list1[len(variance_list1)-1])
 variance_final2.append(variance_list2[len(variance_list2)-1])
 variance_final3.append(variance_list3[len(variance_list3)-1])
# variance_final2.append(1)
# variance_final3.append(1)
 variance_final4.append(1)

 Count_final1.append(Count_list1[len(Count_list1)-1])
 Count_final2.append(Count_list2[len(Count_list2)-1])
 Count_final3.append(Count_list3[len(Count_list3)-1])
#  Count_final4.append(1)
# Count_final2.append(1)
# Count_final3.append(1)
 Count_final4.append(1)

 time_final1.append(time1[len(time1)-1])
 time_final2.append(time2[len(time2)-1])
 time_final3.append(time3[len(time3)-1])
# time_final2.append(1)
# time_final3.append(1)
 time_final4.append(1)










 plotdata.plot(Count_list1, variance_list1,Count_list2, variance_list2,Count_list3, variance_list3,Count_list4, variance_list4, q,L,hz_list,E_time1, time1,E_time2, time2,E_time3, time3,E_time4, time4)
 plotdata.Store(variance_final1, variance_final2,variance_final3,variance_final4, q , L , file, Count_final1, Count_final2, Count_final3,Count_final4, time1, time2, time3,time4,time_final1,time_final2,time_final3,time_final4)
 
#############################################################################################################33

file.close()
















