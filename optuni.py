import numpy as np
import pyUni10 as uni10
import itertools
#import qsampl 
import random
import copy
import time
import mpo
import qsampl
import env
import utils
import optimize

def Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list, Gamma, Max_number_iteratoin_SVD,Max_number_iteratoin_Steepest,Max_number_iteratoin_CG):
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
      L_position=i
      optimize.optimize_function(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list, Count_list,Gamma)
      if i is not L/2 -1:
        Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
#################################################################
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

def variantional_unitary_process(L,L_lay,Model,J,Fieldz,hz_list,Ulist=None,W=8,U_delta=0.,J2=0.,d=2,Randomness='Fixed',Method='SteepestDescentploy',maxiter=50,save_ulist=False):

    utils.print_banner(Model,L,L_lay,J,Fieldz,Method)

    #################INITIALIZATION####################
    Max_SVD_iteratoin = maxiter
    Max_Steepest_iteratoin = maxiter
    Max_CG_iteratoin =  maxiter

    #mpo.intialize_coupling(L, hz_list, W, Randomness, Model) # For Heisenberg
    if len(hz_list) == 0:
        hz_list = qsampl.initialize_hzlist(L, hz_list, Fieldz)
    #qsampl.initialize_hzlist(L, hz_list, Fieldz)
    mpo_list1, mpo_list2, mpo_boundy_list=mpo.make_mpo_H( L, J, hz_list, Fieldz, J2,Model )
    trH2=mpo.contraction_MPO_trH2( mpo_list2, mpo_boundy_list, L )
    #print "Initial energy variance: ", trH2

    if Ulist is None:
        U_list=mpo.intialize_unitary_list( L, L_lay, d, U_delta)
    else:
        U_list = optimize.copy_U_list_function(Ulist, L, L_lay)
    copy_U_list=optimize.copy_U_list_function(U_list, L, L_lay)

    # make a mpo representation of unitary U, up stand for U, down for U^{T} # 
    mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up')
    mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down')
    
    print "****** Making Environment ******"
    perl_label_up=[None]*3
    perl_label_up=[None]*3
    Bond_IN=[None]*3
    List_position=[ 0, 1, (L/2)-1]
    for i in xrange(3):
        L_position=List_position[i]
        perl_label_up[i], Bond_IN[i]= mpo.make_mpo_U_Label( L_position, L_lay, L, 'up')

    Environment_Left=[None]*(L/2)
    Environment_Right=[None]*(L/2)
    Environment_Uni=[None]*(L/2)

    t0 = time.time()
    for i in xrange(1):
        for i in xrange((L/2)-1):
            L_position=i
            Environment_Left[i]=env.Env_left (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
    print "Left Env used   ", time.time() - t0, "seconds wall time"
 
    t0 = time.time()
    for i in xrange(1):
        for i in xrange((L/2)-1):
            L_position=(L/2)-1-i
            Environment_Right[L_position]=env.Env_right (mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Right)
    print "Right Env used  ", time.time() - t0, "seconds wall time"
    
    t0 = time.time()
    for i in xrange(1):
        for i in xrange( (L/2) ):
            L_position=i
            Environment_Uni[L_position]=env.Environment_uni_function(mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left, Environment_Right)
    print "Env Unitary used",time.time() - t0, "seconds wall time"
    
    Env_Uni_inner, Gamma=Initialize_function(L, L_lay)

    for i in xrange(len(L_lay)):
        L_position=0
        L_lay_selected=i
        env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner )
        print (Env_Uni_inner[0][L_lay_selected]*U_list[0][L_lay_selected])[0], 'Left'
    
    for i in xrange(0,len(L_lay),2):
        L_position=(L/2)-1
        L_lay_selected=i
        env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner )
        print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0], 'Right'
    
    if L > 4: 
        for i in xrange(len(L_lay)):
            L_position=1
            L_lay_selected=i
            env.Env_Uni_inner_function(U_list, Environment_Uni, perl_label_up, Bond_IN, L_lay,L_lay_selected, L_position, Env_Uni_inner )
            print (Env_Uni_inner[L_position][L_lay_selected]*U_list[L_position][L_lay_selected])[0],'middle'


    print "****** Starting Optimization ******"
    U_list=copy_U_list
    mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up')
    mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up')
    mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down')
    E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)
    
    E_list4=[]
    Count_list4=[]
    Optimi_full_process(U_list,mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list,Environment_Left,Environment_Right,perl_label_up, Environment_Uni,Env_Uni_inner, Bond_IN,d,L,L_lay,L_position,Method ,Max_SVD_iteratoin, Max_Steepest_iteratoin,Max_CG_iteratoin, E_list4, Count_list4,Gamma,Max_SVD_iteratoin,Max_Steepest_iteratoin,Max_CG_iteratoin)
    
    E_f=mpo.Energy_through_env(U_list, L_lay, L, mpo_U_list_up, mpo_U_list_down, mpo_list2, mpo_boundy_list, L_position, d, Environment_Left)

    variance_Energy=mpo.variance_Energy_function(trH2, E_f, L)

    print "Final energy variance: ", variance_Energy
    
    mpo_U_list_up= mpo.make_mpo_U_list(U_list, L_lay, L, 'up')
    mpo_U_list_down= mpo.make_mpo_U_list(U_list, L_lay, L, 'down')
    
    if save_ulist:
        for i in xrange(len(U_list)):
            for j in xrange(len(U_list[0])):
                U_list[i][j].save("./ulist/U%d%d"%(i,j))

    return U_list
    

