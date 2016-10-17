import pyUni10 as uni10
#import sys
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


L=10             # Numbers of particles
L_lay=[0,1]   #Numbers of layers <= 4
d=2             #pysical bond-dimension  
chi=5           #bond-dimension of MPO
W=8             #random interval, [-W,W]
J=1.0           #coupling
hz_list=[]
U_delta=0.00  #if it's zero, U_list is intialize by Identity
Method='SVD'    #There are three methods: Conjugate, SVD, SteepestDescent
Randomness='Fixed'
Max_number_iteratoin=20   # maximum number of sweeps 
Max_SVD_iteratoin=20     #maximum number of SVD iteration
Max_Steepest_iteratoin=20  # maximum number of SteepestDescent iteration
#######################################################################################3

mpo.intialize_coupling(L, hz_list, W,Randomness)

print hz_list
mpo.print_trH2(L, hz_list,J)

mpo.avarage_Energy_power_2(L, hz_list,J)



######### retrun MPO list, 1: one-site MPO, 2: two-site MPO; mpo_uni10_bl, mpo_uni10_br, mpo_uni10_bl2, mpo_uni10_br2###
mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list )
#print len(mpo_list2)
#mpo.contraction_MPO_trH2( mpo_list1, mpo_list2, mpo_boundy_list )

U_list=mpo.intialize_unitary_list( L/2, L_lay, d, U_delta)
#U_list_Trans=mpo.unitary_list_transpose( U_list, L, L_lay)
#print len(U_list)



#U1=copy.copy(U_list[4][2])
#U1.setLabel([4,5,2,3])
#print U1.transpose()
#print U_list[4][2]*U1.transpose()


L_position=4
############ make a mpo representation of unitary U at position L_position  ###############
example_mpo_U=mpo.make_mpo_U(U_list, L_position, L_lay, L, 'up')

######### make a mpo representation of unitary U, up stand for U, down for U^{T} ######### 
mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up')
mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down')

perl_label_up=[None]*3
perl_label_up=[None]*3
Bond_IN=[None]*3
List_position=[ 0, 1, (L/2)-1]
for i in xrange(3):
 L_position=List_position[i]
 perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up')
#print perl_label_up[0], perl_label_up[1], perl_label_up[2]

L_lay_h=[0]
for i in xrange(len(L_lay)):
 mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay_h, L, 'up')
 mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay_h, L, 'down')
 L_lay_h.append([i+1])
 #mpo_U_list_down[4].setLabel([0,4,2,1])
 #print mpo_U_list_down[0].printDiagram()
 #mpo_U_list_down[8].setLabel([3,1,0])
 #print mpo_U_list_up[8]*mpo_U_list_down[8]

mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up')
mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down')



Environment_Left=[None]*(L/2)
Environment_Right=[None]*(L/2)
Environment_Uni=[None]*(L/2)
t0 = time.time()
for i in xrange((L/2)-1):
 L_position=i
 #print L_position
 Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
 #print Environment_Left[i].printDiagram()
#print time.time() - t0, "seconds wall time"
 
#t0 = time.time()
for i in xrange((L/2)-1):
 L_position=(L/2)-1-i
 #print L_position
 Environment_Right[L_position]=env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Right)
 #print Environment_Right[L_position].printDiagram()
 #print Environment_Right[0]
 #print (Environment_Right[1] * Environment_Left[0])[0] 
#print time.time() - t0, "seconds wall time"
#t0 = time.time()
for i in xrange( (L/2) ):
 L_position=i
 Environment_Uni[L_position]=env.Environment_uni_function(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right)
 #print  (mpo_U_list_up[L_position] * Environment_Uni[L_position])[0]

#print time.time() - t0, "seconds wall time"



Env_Uni_inner=[]
for i in xrange(L/2):
 Env_Uni_inner.append([])
for i in xrange(L/2):
 for j in xrange(len(L_lay)):
  Env_Uni_inner[i].append([None])

for i in xrange(len(L_lay)):
 L_position=0
 L_lay_selected=i
 env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner )

 print (Env_Uni_inner[0][L_lay_selected]*U_list[0][L_lay_selected])[0]

for i in xrange(0,len(L_lay),2):
 L_position=(L/2)-1
 L_lay_selected=i
 env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner )

 print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0]




for i in xrange(len(L_lay)):
 L_position=3
 L_lay_selected=i
 env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner )

 print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0]
 






for j in xrange(Max_number_iteratoin):
 for i in xrange( L/2 ):
  #making Env right
  if i is 0:
   for m in xrange((L/2)-1):
    L_position1=(L/2)-1-m
    #print L_position
    Environment_Right[L_position1]=env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position1, d, Environment_Right)

  print i
  L_position=i
  optimize.optimize_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin)
  if i is not L/2 -1:
   Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)


mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up')
mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up')
mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down')

for i in xrange((L/2)-1):
 L_position=i
 Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
 #print Environment_Left[i].printDiagram()
#print time.time() - t0, "seconds wall time"
 
for i in xrange((L/2)-1):
 L_position=(L/2)-1-i
 #print L_position
 Environment_Right[L_position]=env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Right)

print (Environment_Right[1] * Environment_Left[0])[0] 
print (Environment_Right[2] * Environment_Left[1])[0] 

Avarage_E_power2=(Environment_Right[1] * Environment_Left[0])[0]
variance_Energy=mpo.variance_Energy_function(L, hz_list,J,Avarage_E_power2)
print variance_Energy

