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
################################## Set parameter #######################################
Tech='Regular'                   # MERA, Regular  
Model='Ising'                    #could be: 'Ising' and or 'Heisenberg'
L=8                             #Numbers of particles should be even!
L_lay=[0,1,2,3,4]                    #Numbers of layers <= 5
d=2                              #pysical bond-dimension  
chi=5                            #bond-dimension of MPO
W=8                              #random interval, [-W,W]
J=1.0                            #coupling, 1.0
J2=0.30                           #coupling, 0.3
Fieldz=0.60                       #Field in z direction, 0.6
hz_list=[]                        # list of randomness 
U_delta=0.00                    #if it's zero, U_list is intialized by Identity
Method='SteepestDescent'          #methods: CGarmjo,CGploy, SVD, SteepestDescent, SteepestDescentploy 
Randomness='Fixed'                  #Fixed
Max_number_iteratoin_SVD=30         # maximum number of sweeps for SVD method 
Max_number_iteratoin_Steepest=20   # maximum number of sweeps for SteepestDescent method
Max_number_iteratoin_CG=20         #maximum number of sweeps for CG method
Max_SVD_iteratoin=16               #maximum number of SVD iteration
Max_Steepest_iteratoin=16          #maximum number of SteepestDescent iteration
Max_CG_iteratoin=16                #maximum number of SteepestDescent iteration
#######################################################################################3
def Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list, Gamma,Tech):
 if Method is 'SVD':
  Length=Max_number_iteratoin_SVD
 elif Method is 'SteepestDescent':
  Length=Max_number_iteratoin_Steepest
 elif (Method is 'CGarmjo') or (Method is 'CGpoly') or (Method is 'SteepestDescentploy'):
  Length=Max_number_iteratoin_CG
 for j in xrange(Length):
  for i in xrange( L/2 ):
   #making Env right
   if i is 0:
    for m in xrange((L/2)-1):
     L_position1=(L/2)-1-m
     #print L_position
     Environment_Right[L_position1]=env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position1, d, Environment_Right)

   print i
   L_position=i
   optimize.optimize_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list,Gamma,Tech)
   if i is not L/2 -1:
    Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
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
mpo.intialize_coupling(L, hz_list, W,Randomness, Model)

print hz_list

##########    only works for Heisenberg    ##########################
#trH2=mpo.print_trH2(L, hz_list,J)
#print trH2 
#mpo.avarage_Energy_power_2(L, hz_list,J)

######### retrun MPO list, 1: one-site MPO, 2: two-site MPO; mpo_uni10_bl, mpo_uni10_br, mpo_uni10_bl2, mpo_uni10_br2###
mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model )
trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L )
print trH2

U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta)
#U_list_Trans=mpo.unitary_list_transpose( U_list, L, L_lay)
#print len(U_list)
copy_U_list=optimize.copy_U_list_function(U_list, L, L_lay)
copy_U_list1=optimize.copy_U_list_function(U_list, L, L_lay)
copy_U_list2=optimize.copy_U_list_function(U_list, L, L_lay)
copy_U_list3=optimize.copy_U_list_function(U_list, L, L_lay)

######### make a mpo representation of unitary U, up stand for U, down for U^{T} ######### 
#print U_list[0][0]
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


#t0 = time.time()
for i in xrange(1):
 for i in xrange((L/2)-1):
  L_position=i
  Environment_Left[i]=env.Env_left(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
  #print Environment_Left[i].printDiagram()
#print Environment_Left[0].printDiagram()
#print time.time() - t0, "seconds wall time"
 
 


t0 = time.time()
for i in xrange(1):
 for i in xrange((L/2)-1):
  L_position=(L/2)-1-i
  Environment_Right[L_position]=env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Right)
  #print Environment_Right[L_position].printDiagram()
 
#print Environment_Right[0]
print (Environment_Right[1] * Environment_Left[0])[0] 
print time.time() - t0, "seconds wall time"




t0 = time.time()
for i in xrange(1):
 for i in xrange( (L/2) ):
  L_position=i
  Environment_Uni[L_position]=env.Environment_uni_function(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right)
  print  (mpo_U_list_up[L_position] * Environment_Uni[L_position])[0]

print time.time() - t0, "seconds wall time"


Env_Uni_inner, Gamma=Initialize_function(L, L_lay)


if Tech is 'Regular':
 for i in xrange(len(L_lay)):
  L_position=0
  L_lay_selected=i
  env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner,Tech )

  print (Env_Uni_inner[0][L_lay_selected]*U_list[0][L_lay_selected])[0], 'Left'
elif Tech is 'MERA':
 for i in xrange(len(L_lay)-3):
  L_position=0
  L_lay_selected=i
  env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner,Tech )

  print (Env_Uni_inner[0][L_lay_selected]*U_list[0][L_lay_selected])[0], 'Left'


if Tech is 'Regular':
 for i in xrange(0,len(L_lay),2):
  L_position=(L/2)-1
  L_lay_selected=i
  env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner,Tech )

  print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0], 'Right'
elif Tech is 'MERA':
 for i in xrange(0,len(L_lay)-3,2):
   L_position=(L/2)-1
   L_lay_selected=i
   env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner,Tech )

   print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0], 'Right'

if L > 4: 
 for i in xrange(len(L_lay)-5):
  L_position=1
  L_lay_selected=i
  env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner,Tech )

  print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0],'middle'
 
 for i in xrange(len(L_lay)):
  L_position=2
  L_lay_selected=i
  env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner,Tech )

  print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0],'middle'


##################################### Stepest Method ###################################
#Method='SteepestDescent'
#print 'hi',Method
#E_list=[]
#Count_list=[]
#Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list,Gamma )


##E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Tech)
##print  E_f

#print Count_list
####################################### SVD ######################################


Method='SVD'
U_list=copy_U_list
mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down', Tech)
E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
print  E_f

E_list1=[]
Count_list1=[]
Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list1, Count_list1,Gamma,Tech )




#print Count_list1
#E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
#print  E_f

#variance_Energy=mpo.variance_Energy_function(trH2, E_f, L)
#print variance_Energy
############################################   CG-armjo   ###################################################


Tech='MERA'
L_lay=[0,1,2,3,4,5,6,7]
perl_label_up=[None]*4
perl_label_up=[None]*4
Bond_IN=[None]*4
List_position=[ 0, 1 , 2, (L/2)-1]
for i in xrange(4):
 L_position=List_position[i]
 perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up',  Tech)
#print perl_label_up[0], perl_label_up[1], perl_label_up[2]
Env_Uni_inner, Gamma=Initialize_function(L, L_lay)
U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta)
copy_U_list1=optimize.copy_U_list_function(U_list, L, L_lay)
U_list=copy_U_list1


mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down', Tech)
E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
print  E_f
print  E_list1[0], 'must be equall'



Method='SVD'
E_list2=[]
Count_list2=[]
Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list2, Count_list2,Gamma,Tech )
#print Count_list2
E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
print  E_f




#variance_Energy=mpo.variance_Energy_function(trH2, E_f, L)
#print variance_Energy
##############################################################################################################33

#############################################   CG-armjo   #########################################
#print 'CGpoly','\n','\n','\n'
#U_list=copy_U_list2
#mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
#mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
#mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down', Tech)
#E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Tech)
#print  E_f
#print  E_list[0], 'must be equall'

#Method='CGpoly'
#E_list3=[]
#Count_list3=[]
#Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list3, Count_list3,Gamma,Tech )

#print Count_list3

#E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left,Tech)
#print  E_f

#variance_Energy=mpo.variance_Energy_function(trH2, E_f, L)
#print variance_Energy

###########################################################
#print 'SteepestDescentploy','\n','\n','\n'
#U_list=copy_U_list3
#mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
#mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up', Tech)
#mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down', Tech)
#E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Tech)
#print  E_f
#print  E_list[0], 'must be equall'


#Method='SteepestDescentploy'
#E_list4=[]
#Count_list4=[]
#Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list4, Count_list4,Gamma )

#print Count_list4

#E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Tech)
#print  E_f

#variance_Energy=mpo.variance_Energy_function(trH2, E_f, L)
#print variance_Energy








#############################################################################################################33



#variance_list=[]
#variance_list=[ (trH2-E_list[i])*(1.0/(2**L)) for i in xrange(len(E_list)) ]
variance_list1=[]
variance_list1=[ (trH2-E_list1[i])*(1.0/(2**L)) for i in xrange(len(E_list1)) ]
variance_list2=[]
variance_list2=[ (trH2-E_list2[i])*(1.0/(2**L)) for i in xrange(len(E_list2)) ]
#variance_list3=[]
#variance_list3=[ (trH2-E_list3[i])*(1.0/(2**L)) for i in xrange(len(E_list3)) ]
#variance_list4=[]
#variance_list4=[ (trH2-E_list4[i])*(1.0/(2**L)) for i in xrange(len(E_list4)) ]

print variance_list1[len(variance_list1)-1],variance_list2[len(variance_list2)-1]


#plt.plot(Count_list,variance_list, 'ro',label='SteepestDescent-armjo')
plt.plot( Count_list1,variance_list1,'bs',label='Regular-architecture' )
plt.plot( Count_list2, variance_list2,'g.',label='MERA' )
#plt.plot( Count_list3, variance_list3,'yh',label='CG-poly' )
#plt.plot( Count_list4, variance_list4,'c+',label='SteepestDescent-poly' )

plt.xlabel('Iterations', fontsize=20)
plt.ylabel('$\sigma$', fontsize=20)
plt.legend(loc='upper right')
plt.savefig('Figs/Convergance.pdf')
plt.show()
#print Gamma


#a=[0, 1,2]
#b=[-4,5,6]

#file = open("list.txt", "w")
#for index in range(len(a)):
#    file.write(str(a[index]) + " " + str(b[index]) + "\n")
#file.close()
#R=loadtxt('list.txt')


