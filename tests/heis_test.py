import pyUni10 as uni10
#import syss
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import numpy as np
#matplotlib.use('Agg')
import sys
sys.path.append('/Users/sunchong/CodingDir/TNMBL')
import qsampl 
import utils
import optuni
import random
import copy
import time
import mpo
import env
import optimize
#np.set_printoptions(linewidth=1000)
################################## Set parameter #######################################
Model='Heisenberg'               #could be: 'Ising' and or 'Heisenberg'
L=int(sys.argv[1])               #Numbers of particles should be even!
llayer = int(sys.argv[2])
read_u = False
if llayer == 1:
    L_lay = [0]
elif llayer == 2:
    L_lay = [0,1]
else:
    raise Exception("Layer number can only be 1 or 2!")

d=2                              #pysical bond-dimension  
chi=5                            #bond-dimension of MPO

J=1                           #coupling, 1.0
Fieldz=1.0                            #Field in z direction, 0.6


W=8                              #random interval, [-W,W]
J2=0.0                           #coupling, 0.3
hz_list=[]                        # list of randomness 
U_delta=0.001                   #if it's zero, U_list is intialize by Identit
Method='CGploy'          #methods: CGarmjo,CGploy, SVD, SteepestDescent, SteepestDescentploy 
randomhz = True
Randomness='Fixed'
Max_number_iteratoin_SVD=60        # maximum number of sweeps for SVD method 
Max_number_iteratoin_Steepest=60   # maximum number of sweeps for SteepestDescent method
Max_number_iteratoin_CG=60   # maximum number of sweeps for SteepestDescent method
Max_SVD_iteratoin=30               #maximum number of SVD iteration
Max_Steepest_iteratoin=30          #maximum number of SteepestDescent iteration
Max_CG_iteratoin=30          #maximum number of SteepestDescent iteration
######################################################


#U_list = optuni.variantional_unitary_process(L,L_lay,Model,J,Fieldz,hz_list,U_delta=U_delta)

hz_list = qsampl.initialize_hzlist(L,hz_list,Fieldz,W,randomhz,Randomness)

if read_u:
  try:
    U_list = utils.load_ulist(L/2,llayer,outdir='./ulist/'+Model+'/J%0.1fh%0.1f/'%(J,Fieldz))
  except:
    U_list = optuni.variantional_unitary_process(L,L_lay,Model,J,Fieldz,hz_list,U_delta=U_delta,randomhz=randomhz,Randomness=Randomness)
    utils.save_ulist(U_list,outdir='./ulist/'+Model+'/J%0.1fh%0.1f/'%(J, Fieldz))
else:
  U_list = optuni.variantional_unitary_process(L,L_lay,Model,J,Fieldz,hz_list,U_delta=U_delta,randomhz=randomhz,Randomness=Randomness)
  utils.save_ulist(U_list,outdir='./ulist/'+Model+'/J%0.1fh%0.1f/'%(J, Fieldz))


  


## Iterative
#U_list0 = optuni.variantional_unitary_process(L,L_lay,Model,0.1,1,hz_list,U_delta=0.01)
#U_list1 = optuni.variantional_unitary_process(L,L_lay,Model,0.4,1,hz_list,Ulist=U_list0,U_delta=0.01)
#U_list2 = optuni.variantional_unitary_process(L,L_lay,Model,0.8,1,hz_list,Ulist=U_list1,U_delta=0.01)
#U_list = optuni.variantional_unitary_process(L,L_lay,Model,1,1,hz_list,Ulist=U_list2,U_delta=0.1)
#
#
#print U_list[0][0]
#print U_list3[0][0].getBlock()
#
#exit()
#
#

def spectrum():
    Egs, v = qsampl.heisenberg_ED(L,J,hz_list) 
    np.savetxt("./data/spectrum/%s_L%d_ed_spectrum.txt"%(Model,L), Egs)
    spect = qsampl.get_spectrum(U_list,Model,L,hzlist=hz_list,outdir="./data/spectrum/")  
 
   
def obtain_data():

    Egs, v = qsampl.heisenberg_ED(L,J,hz_list)
    psi = utils.vec2uni(v,L)
    psi = qsampl.UdagPsi(U_list,psi,L)
    ave,err = qsampl.measure_pstr(U_list,L,psi,Model,J,Fieldz,hz_list,analys_coefs=True,coutdir='./data/coefs/',eoutdir='./data/errors/',tol=5e-1)
    print "True ED: ", Egs[0]
    print "ED from Paulis: ", ave
    print "qsample err: ", err
    
    oerr = qsampl.heisenberg_variance(L,J,hz_list)
    print "original error: ", oerr

   

if __name__ == "__main__":
    obtain_data()
    spectrum()    




