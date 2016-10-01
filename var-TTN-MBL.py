import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random 
import copy
import MPO

L=8
L_lay=[0,1,2]  #three layers
d=2
chi=5
J=5
hz_list=[]
for i in xrange(L):
 hz_list.append(np.random.uniform(-J,J))

print hz_list
hz_list_2=[x**2 for x in hz_list]
print hz_list_2
print sum(hz_list_2)
trH2=(2.00**L)*(3.0/16.00*(L-1)+1.0/4.0*sum(hz_list_2))
print  trH2

mpo_list1, mpo_list2, mpo_bl1 , mpo_br1, mpo_bl2, mpo_br2=MPO.make_mpo_H(L, J, hz_list )
MPO.contraction_MPO_trH2(mpo_list1, mpo_list2, mpo_bl1 , mpo_br1, mpo_bl2, mpo_br2 )
