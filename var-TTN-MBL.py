import pyUni10 as uni10
#import sys
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random 
import copy
import MPO

L=14
L_lay=[0,1,2,3,4]  #Numbers of layers <= 4
d=2
chi=5
W=5
J=3.0
hz_list=[]
U_delta=0.400  #if it's zero, U_list is intialize by Identity



MPO.intialize_coupling(L, hz_list, W)
MPO.print_trH2(L, hz_list,J)




# retrun MPO list, 
mpo_list1, mpo_list2, mpo_boundy_list=MPO.make_mpo_H( L, J, hz_list )

print len(mpo_list2)


MPO.contraction_MPO_trH2( mpo_list1, mpo_list2, mpo_boundy_list )

U_list=MPO.intialize_unitary_list( L, L_lay, d, U_delta)
U_list_Trans=MPO.unitary_list_transpose( U_list, L, L_lay)




#U1=copy.copy(U_list[4][2])
#U1.setLabel([4,5,2,3])
#print U1.transpose()
#print U_list[4][2]*U1.transpose()
 
L_position=6

example_mpo_U=MPO.make_mpo_U(U_list,L_position, L_lay, L, 'up')

mpo_U_list_up= MPO.make_mpo_U_list(U_list, L_lay, L, 'up')
mpo_U_list_down= MPO.make_mpo_U_list(U_list, L_lay, L, 'down')


example_mpo_Up=MPO.make_mpo_U(U_list,L_position, L_lay, L, 'down')

L_lay_h=[0]
for i in xrange(5):
 mpo_U_list_up= MPO.make_mpo_U_list(U_list, L_lay_h, L, 'up')
 mpo_U_list_down= MPO.make_mpo_U_list(U_list, L_lay_h, L, 'down')
 L_lay_h.append([i+1])
 mpo_U_list_down[5].setLabel([0,4,2,1])
 #print mpo_U_list_down[0].printDiagram()
 #mpo_U_list_down[5].setLabel([3,1,0])


 print mpo_U_list_up[5]*mpo_U_list_down[5]






