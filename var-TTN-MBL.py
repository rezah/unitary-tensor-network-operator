import pyUni10 as uni10
#import sys
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random 
import copy
import MPO
import Env
L=20
L_lay=[0,1,2,3]  #Numbers of layers <= 4
d=2
chi=5
W=5
J=3.0
hz_list=[]
U_delta=0.00  #if it's zero, U_list is intialize by Identity
#######################################################################################3


MPO.intialize_coupling(L, hz_list, W)
MPO.print_trH2(L, hz_list,J)

MPO.variance_Energy(L, hz_list,J)



######### retrun MPO list, 1: one-site MPO, 2: two-site MPO; mpo_uni10_bl, mpo_uni10_br, mpo_uni10_bl2, mpo_uni10_br2###
mpo_list1, mpo_list2, mpo_boundy_list=MPO.make_mpo_H( L, J, hz_list )

#MPO.contraction_MPO_trH2( mpo_list1, mpo_list2, mpo_boundy_list )

U_list=MPO.intialize_unitary_list( L, L_lay, d, U_delta)
#U_list_Trans=MPO.unitary_list_transpose( U_list, L, L_lay)




#U1=copy.copy(U_list[4][2])
#U1.setLabel([4,5,2,3])
#print U1.transpose()
#print U_list[4][2]*U1.transpose()


L_position=4
############ make a mpo representation of unitary U at position L_position  ###############
example_mpo_U=MPO.make_mpo_U(U_list, L_position, L_lay, L, 'up')

######### make a mpo representation of unitary U, up stand for U, down for U^{T} ######### 
mpo_U_list_up= MPO.make_mpo_U_list(U_list, L_lay, L, 'up')
mpo_U_list_down= MPO.make_mpo_U_list(U_list, L_lay, L, 'down')

L_lay_example=[0,1,2]
L_position=0
perl_label_up, Bond_IN= MPO.make_mpo_U_Label( L_position, L_lay_example, L, 'up')
#print perl_label_up, Bond_IN

L_lay_h=[0]
for i in xrange(3):
 mpo_U_list_up= MPO.make_mpo_U_list(U_list, L_lay_h, L, 'up')
 mpo_U_list_down= MPO.make_mpo_U_list(U_list, L_lay_h, L, 'down')
 L_lay_h.append([i+1])
 #mpo_U_list_down[4].setLabel([0,4,2,1])
 #print mpo_U_list_down[0].printDiagram()
 #mpo_U_list_down[8].setLabel([3,1,0])
 #print mpo_U_list_up[8]*mpo_U_list_down[8]

mpo_U_list_up= MPO.make_mpo_U_list(U_list, L_lay, L, 'up')
mpo_U_list_down= MPO.make_mpo_U_list(U_list, L_lay, L, 'down')



Environment_Left=[None]*(L/2)
Environment_Right=[None]*(L/2)
Environment_Uni=[None]*(L/2)

for i in xrange((L/2)-1):
 L_position=i
 Environment_Left[i]=Env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
 #print Environment_Left[i].printDiagram()
 
 

for i in xrange((L/2)-1):
 L_position=(L/2)-1-i
 #print L_position
 Environment_Right[L_position]=Env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Right)
 #print Environment_Right[L_position].printDiagram()
 #print Environment_Right[0]
print (Environment_Right[1] * Environment_Left[0])[0] 
print (Environment_Right[8] * Environment_Left[7])[0]
print (Environment_Right[3] * Environment_Left[2])[0]

L_position=0
Environment_Uni[L_position]=Env.Environment_uni(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right)





