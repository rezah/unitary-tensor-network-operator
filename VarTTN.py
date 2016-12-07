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
Tech='MERA'                  # MERA, Regular  
Model='Heisenberg'              #could be: 'Ising' and or 'Heisenberg'
L=8+6*(2*1)                            #Numbers of particles should be 4*k, k={1,2,3,4,...}!6*(2k)+8
L_lay=[0,1,2,3,4,5,6,7,8]               #Numbers of layers <= 5
realizations=1                  #Number of realizations        
Accuracy=8.00e-5                #Accuracy of variational variance  
d=2                              #pysical bond-dimension  
chi=5                            #bond-dimension of MPO
W=8                              #random interval, [-W,W]
J=1.0                            #coupling, 1.0
J2=0.30                          #coupling, 0.3
Fieldz=0.60                      #Field in z direction, 0.6
hz_list=[]                       # list of randomness 
U_delta=0.00                     #if it's zero, U_list is intialized by Identity
Method='SteepestDescent'         #methods: CGarmjo, CGpoly, SVD, SteepestDescent, SteepestDescentploy 
Randomness='Fixed'               #Fixed
Max_number_iteratoin_SVD=1       # maximum number of sweeps for SVD method 
Max_number_iteratoin_Steepest=400  # maximum number of sweeps for SteepestDescent method
Max_number_iteratoin_CG=400         #maximum number of sweeps for CG method

Max_SVD_iteratoin=1               #maximum number of SVD iteration
Max_Steepest_iteratoin=20          #maximum number of SteepestDescent iteration
Max_CG_iteratoin=20                #maximum number of SteepestDescent iteration
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
  print 'count=', count, 'Accuracy=', abs((E0-E1)/E1),'E0=', E0,'E1=', E1, '  ', trH2, 'Loop_iter=', Loop_iter, 'Count_list=', Count_list[len(Count_list)-1], 'Variance=', (trH2-E1) * (1.0/(2**(L_fixed)))
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
trH2=mpo.print_trH2(L, hz_list,J)
print 'trH2:=', trH2 
variance_Energy=mpo.avarage_Energy_power_2(L, hz_list,J)
print variance_Energy
######### retrun MPO list, 1: one-site MPO, 2: two-site MPO; mpo_uni10_bl, mpo_uni10_br, mpo_uni10_bl2, mpo_uni10_br2###########################################
mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2, Model, Tech, L_lay )

#print mpo_list2[0].printDiagram(),mpo_list2[1].printDiagram(),mpo_list2[2].printDiagram(),mpo_list2[3].printDiagram(),mpo_list2[4].printDiagram(),mpo_list2[5].printDiagram()

trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L, Tech, L_lay )
print 'trH2:=', trH2

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

t0=time.time()
for i in xrange((L_m/2)-1):
 L_position=i
 print 'i',i
 Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
print time.time() - t0, "Seconds, Left"

t0=time.time()
for i in xrange((L_m/2)-1):
 L_position=(L_m/2)-1-i
 print L_position
 Environment_Right[L_position]=env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Right)

print time.time() - t0, "Seconds, Right"

print (Environment_Right[1] * Environment_Left[0])[0] 

t0=time.time()
for i in xrange(1):
 for j in xrange(2):
  L_position=0
  print 'i', i
  L_lay_selected=j
  Environment_Uni[L_position]=env.Environment_uni_function(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right) 
  env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner,Tech )
  print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0], 'L'
print time.time() - t0, "Seconds, L"


t0=time.time()
for i in xrange(1):
 for j in xrange(8):
  L_position=1
  print 'i=', i, j
  L_lay_selected=j
  Environment_Uni[L_position]=env.Environment_uni_function(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right) 
  env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner,Tech )
  print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0], 'M'
print time.time() - t0, "Seconds, M"


t0=time.time()
for i in xrange(1):
 for j in xrange(3):
  L_position=3
  print 'j=', i
  L_lay_selected=j
  Environment_Uni[L_position]=env.Environment_uni_function(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right) 
  env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner,Tech )
  print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0], 'R'
print time.time() - t0, "Seconds, R"

#Environment_Left=[None]*(L_m/2)
#Environment_Right=[None]*(L_m/2)
#Environment_Uni=[None]*(L_m/2)
#Env_Uni_inner, Gamma=Initialize_function(L_m, L_lay)



file = open("Data/varianceAll.txt", "w")
variance_final1=[]
variance_final2=[]
variance_final3=[]
variance_final4=[]

Count_final1=[]
Count_final2=[]
Count_final3=[]
Count_final4=[]




for q in xrange(realizations):
 print '\n', 'q_MERA_ternary=', q, '\n'
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
 
 
 Randomness='Fixed'
 mpo.intialize_coupling(L, hz_list, W, Randomness, Model)
 
 trH2=mpo.print_trH2(L, hz_list,J)
 print 'trH2=', trH2 
 
 
 
 
 
 Tech='MERA'
 L_lay=[0,1,2,3,4,5,6,7,8]
 mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model,Tech,L_lay )
 trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L, Tech,L_lay)
 print '\n', 'trH2=', trH2
 Method='SVD'
 U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta,Tech)
 mpo_U_list_up = mpo.make_mpo_U_list(U_list, L_lay, (((L-8)/3)+4), 'up', Tech)
 mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, (((L-8)/3)+4), 'down', Tech)



 List_position=[ 0, 1, 2, ((((L-8)/3)+4)/2)-1]
 for i in xrange(4):
  L_position=List_position[i]
  perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)

 Environment_Left=[None]*((((L-8)/3)+4)/2)
 Environment_Right=[None]*((((L-8)/3)+4)/2)
 Environment_Uni=[None]*((((L-8)/3)+4)/2)
 Env_Uni_inner, Gamma=Initialize_function((((L-8)/3)+4), L_lay)
 E_f=mpo.Energy_through_env(U_list, L_lay, (((L-8)/3)+4), mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
 print 'E_0 = ', E_f,'\n'
 
 Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,(((L-8)/3)+4),L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list4, Count_list4,Gamma,Tech,trH2,E_time4, time4)
 

 
 
#####################################################################################
 
 
 
 Tech='MERA'
 L_lay=[0,1,2,3,4,5]
 mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model,Tech,L_lay )
 trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L, Tech,L_lay)
 print '\n', 'trH2=', trH2
 Method='SVD'
 U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta,Tech)
 mpo_U_list_up = mpo.make_mpo_U_list(U_list, L_lay, L/2, 'up', Tech)
 mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L/2, 'down', Tech)



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
 
 Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L/2,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list3, Count_list3,Gamma,Tech,trH2,E_time1, time1)
 

 #######################################################################################
 Tech='Regular'
 L_lay=[0,1,2,3,4]
 U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta,Tech)
 mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model,Tech,L_lay )
 trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L, Tech,L_lay)
 print '\n', 'trH2=', trH2

 List_position=[ 0, 1, 2, (L/2)-1]
 for i in xrange(4):
  L_position=List_position[i]
  perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)

 Environment_Left=[None]*(L/2)
 Environment_Right=[None]*(L/2)
 Environment_Uni=[None]*(L/2)
 Env_Uni_inner, Gamma=Initialize_function(L, L_lay)
 E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
 print ' E_0= ',  E_f, '\n'


 Method='SVD'
 U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta,Tech)
 mpo_U_list_up = mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
 mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down', Tech)
 Env_Uni_inner, Gamma=Initialize_function(L, L_lay)
 for i in xrange(4):
  L_position=List_position[i]
  perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)
 E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
 Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list1, Count_list1,Gamma,Tech,trH2,E_time2, time2 )
 

 ###################################   MERA  ##########################
 print '\n', 'q_MERA=', q, '\n'
 Tech='MERA'
 L_lay=[0,1,2,3,4,5,6,7]
 for i in xrange(4):
  L_position=List_position[i]
  perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)
 Env_Uni_inner, Gamma=Initialize_function(L, L_lay)
 U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta,Tech)
 mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model,Tech,L_lay )
 trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L, Tech,L_lay)
 print '\n', 'trH2=', trH2
 mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
 mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down', Tech)
 E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
 print ' E_0= ', E_f,'\n'


 Method='SVD'

 Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list2, Count_list2,Gamma,Tech,trH2,E_time3, time3 )

 
 variance_list1=[ (trH2-E_list1[i])*(1.0/(2**L)) for i in xrange(len(E_list1)) ]
 variance_list2=[ (trH2-E_list2[i])*(1.0/(2**L)) for i in xrange(len(E_list2)) ]
 variance_list3=[ (trH2-E_list3[i])*(1.0/(2**L)) for i in xrange(len(E_list3)) ]
 variance_list4=[ (trH2-E_list4[i])*(1.0/(2**L)) for i in xrange(len(E_list4)) ]

 variance_final1.append(variance_list1[len(variance_list1)-1])
 variance_final2.append(variance_list2[len(variance_list2)-1])
 variance_final3.append(variance_list3[len(variance_list3)-1])
 variance_final4.append(variance_list4[len(variance_list4)-1])

 Count_final1.append(Count_list1[len(Count_list1)-1])
 Count_final2.append(Count_list2[len(Count_list2)-1])
 Count_final3.append(Count_list3[len(Count_list3)-1])
 Count_final4.append(Count_list4[len(Count_list4)-1])

 plotdata.plot(Count_list1, variance_list1,Count_list2, variance_list2,Count_list3, variance_list3,Count_list4, variance_list4, q,L,hz_list,E_time1, time1,E_time2, time2,E_time3, time3,E_time4, time4)
 plotdata.Store(variance_final1, variance_final2,variance_final3,variance_final4, q , L , file, Count_final1, Count_final2, Count_final3,Count_final4, time1, time2, time3,time4)
 
#############################################################################################################33

file.close()
















