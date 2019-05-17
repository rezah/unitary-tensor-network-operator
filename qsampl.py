import pyUni10 as uni10
import copy
import itertools
import numpy as np
import pickle
#np.set_printoptions(linewidth=1000)

def matSx():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [0.0, 1.0, 1.00, 0.0])
  return Mat 


def matSz():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [1.0, 0, 0, -1.0]);
  return Mat 

def matSy():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [0.0, -1.00, 1.00, 0.00]);
  return Mat 


def matIden():
    spin_t=0.5
    dimT = int(2*spin_t + 1)
    Mat=uni10.Matrix(dimT, dimT,[1,0,0,1])
    return Mat

####################UHU####################UHU####################UHU
def get_Hbond_tensor(model,bond_idx,L,J=1.0,Fieldz=1.0,hzlist=[],bond_dim=2,\
                     avebond=True):

    #TODO put Z terms only into the even bonds

    '''
    _________________________________________________________________
    Ising Model:                |    Heisenberg Model 
    -J XX - h Z    (h=Fieldz)   |    J(XX + ZZ - YY) - h Z (h=hzlist)
    ____________________________|____________________________________
                        Result in a tensor:
                              _____
                         0 __|     |__ 2
                         1 __|  H  |__ 3
                             |_____| 
    _________________________________________________________________

    '''

    # make matrix
    sx   = matSx()
    sy   = matSy()
    sz   = matSz()
    s0   = matIden()
    if (model == 'Ising'):
        Hbond_mat  = (-1.*J)*uni10.otimes(sx,sx)
        if avebond:
            if (bond_idx == 0):
                Hbond_mat += (-1.0*Fieldz)*uni10.otimes(sz,s0) 
                Hbond_mat += (-0.5*Fieldz)*uni10.otimes(s0,sz)
            elif (bond_idx == L-2):
                Hbond_mat += (-0.5*Fieldz)*uni10.otimes(sz,s0) 
                Hbond_mat += (-1.0*Fieldz)*uni10.otimes(s0,sz)
            else:
                Hbond_mat += (-0.5*Fieldz)*uni10.otimes(sz,s0) 
                Hbond_mat += (-0.5*Fieldz)*uni10.otimes(s0,sz)

        else:
            if (bond_idx%2 == 0):
                Hbond_mat += (-1.0*Fieldz)*uni10.otimes(sz,s0) 
                Hbond_mat += (-1.0*Fieldz)*uni10.otimes(s0,sz)

        
    elif (model == 'Heisenberg'):
        #TODO try evenly distributed bond
        Hbond_mat  =  (1.*J)*uni10.otimes(sx,sx)
        Hbond_mat += (-1.*J)*uni10.otimes(sy,sy)
        Hbond_mat +=  (1.*J)*uni10.otimes(sz,sz)
        if avebond:
            if (bond_idx == 0):
                Hbond_mat += (-1. *hzlist[0])*uni10.otimes(sz,s0)
                Hbond_mat += (-0.5*hzlist[1])*uni10.otimes(s0,sz)
            elif (bond_idx == L-2):
                Hbond_mat += (-0.5*hzlist[L-2])*uni10.otimes(sz,s0)
                Hbond_mat += (-1. *hzlist[L-1])*uni10.otimes(s0,sz)
            else:
                Hbond_mat += (-0.5*hzlist[bond_idx])*uni10.otimes(sz,s0)
                Hbond_mat += (-0.5*hzlist[bond_idx+1])*uni10.otimes(s0,sz)
        else:
            if (bond_idx%2 == 0):
                Hbond_mat += (-1.*hzlist[bond_idx])*uni10.otimes(sz,s0)
                Hbond_mat += (-1.*hzlist[bond_idx+1])*uni10.otimes(s0,sz)
    else:
        raise Exception("The model can only be Ising or Heisenberg!")


    # construct tensor
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)

    Hbond = uni10.UniTensor([bdi,bdi,bdo,bdo],'Hbond')
    Hbond.putBlock(Hbond_mat)
    Hbond.setLabel([0,1,2,3])

    
    return Hbond


def get_umpo_up(Ulist,L,bond_idx,bond_dim=2):

    l_layer = len(Ulist[0])
    if (l_layer == 1): #will not be called
        U = get_umpo_up_2layers(Ulist,L,bond_idx,bond_dim)
    
    elif (l_layer == 2):
        U = get_umpo_up_2layers(Ulist,L,bond_idx,bond_dim)

    return U

def get_umpo_up_1layer(Ulist,L,bond_idx,bond_dim=2):
    # Will not be called
    if (bond_idx%2 == 0):
        U = Ulist[bond_idx/2]
    else:
        bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
        bdo = uni10.Bond(uni10.BD_OUT, bond_dim)
        U = uni10.UniTensor([bdi,bdi,bdo,bdo],'U0')
        I = matIden()
        U.putBlock(I)
    return U
    

def get_umpo_up_2layers(Ulist,L,bond_idx,bond_dim=2):

    #TODO Combine some tensors in advance and save them to speed up
    # Ulist is a L/2 x 2 list, 2 is the number of layers 
    # bond_idx ranges from 0 to L-2
    l_layer = len(Ulist[0])
    if l_layer != 2:
        raise ValueError("Please make sure your unitary list has TWO layers!")
    NU0 = L/2

    if (bond_idx == 0):   # Left Edge
        '''
              _____
         2 __|     |__ 5
         3 __|  U  |__ 6
         4 __|     |__ 7 
             |_____| 
    
        Where 2,3 will combine with the Hamiltonian,
        4 will combine with U_down
        '''
        U0 = copy.copy(Ulist[0][0])
        U1 = copy.copy(Ulist[0][1])
        U0.setLabel([0,1,2,3])
        U1.setLabel([3,4,5,6])
        U = U0*U1
        U.permute([0,1,4,2,5,6],3)
        U.setLabel([2,3,4,5,6,7])

    elif (bond_idx == L-2):   # Right Edge
        '''
              _____
         4 __|     |__ 5
         2 __|  U  |__ 6
         3 __|     |__ 7 
             |_____| 
       
        Where 2,3 will combine with the Hamiltonian,
        4 will combine with U_down
        '''
        U0 = copy.copy(Ulist[NU0-1][0])
        U1 = copy.copy(Ulist[NU0-2][1])
        U0.setLabel([0,1,2,3])
        U1.setLabel([4,2,5,6])
        U = U0*U1
        U.permute([4,0,1,5,6,3],3)
        U.setLabel([4,2,3,5,6,7])
    else:

        '''
            _____
       4 __|     |__ 5
       2 __|     |__ 6
       3 __|  U  |__ 7 
      20 __|     |__ 8        
           |_____| 
     
        Where 2,3 will combine with the Hamiltonian,
        4,20 will combine with U_down
        '''
        if(bond_idx%2 == 0):   # Even Bond
            l = bond_idx/2
            U0  = copy.copy(Ulist[l][0])
            U1L = copy.copy(Ulist[l-1][1])
            U1R = copy.copy(Ulist[l][1])
            U0.setLabel([0,1,2,3])
            U1L.setLabel([4,2,5,6])
            U1R.setLabel([3,7,8,9])
            U = (U1L*U0)*U1R
            U.permute([4,0,1,7,5,6,8,9],4)
            U.setLabel([4,2,3,20,5,6,7,8]) 

        else:                # Odd Bond
            l = bond_idx//2
            U0L = copy.copy(Ulist[l][0])
            U0R = copy.copy(Ulist[l+1][0])
            U1  = copy.copy(Ulist[l][1])
            U0L.setLabel([0,1,2,3])
            U0R.setLabel([4,5,6,7])
            U1.setLabel([3,6,8,9])
            U = U0L*(U1*U0R)
            U.permute([0,1,4,5,2,8,9,7],4)
            U.setLabel([4,2,3,20,5,6,7,8]) 
    
    return U

def contract_UHU(L, Ulist, model='Ising', J=1.0, Fieldz=1.0, hzlist=[]):
    '''
    Return a list of U^+HU.
    L is the lattice length.
    '''
    #TODO Consider when U has imaginary parts
    #TODO take construction of Hbond out of the loop for ising model

    l_layer = len(Ulist[0])

    UHU = []
    if (l_layer == 1):
      for bidx in xrange(L-1):
        Hbond = get_Hbond_tensor(model,bidx,L,J,Fieldz,hzlist)
        Hbond.setLabel([0,1,2,3])
        if (bidx%2==0):
          uidx = bidx/2
          U_up = copy.copy(Ulist[uidx][0])
          U_dn = copy.copy(U_up)
          U_dn.permute([2,3,0,1],2)
          U_up.setLabel([2,3,4,5])
          U_dn.setLabel([6,7,0,1])
          Hbond_n = (U_dn*Hbond)*U_up
          Hbond_n.permute([6,7,4,5],2)
          Hbond_n.setLabel([0,1,2,3])
          UHU.append(Hbond_n)
        else:
          uidx = bidx//2
          U_up_l = copy.copy(Ulist[uidx][0])
          U_up_r = copy.copy(Ulist[uidx+1][0])
          U_dn_l = copy.copy(U_up_l)
          U_dn_r = copy.copy(U_up_r)
          U_dn_l.permute([2,3,0,1],2)
          U_dn_r.permute([2,3,0,1],2)
          U_dn_l.setLabel([4,5,6,0])
          U_dn_r.setLabel([7,8,1,9])
          U_up_l.setLabel([6,2,10,11])
          U_up_r.setLabel([3,9,12,13])
          Hbond_n = U_dn_l*((U_dn_r*Hbond)*U_up_l)*U_up_r
          Hbond_n.setLabel([0,1,2,3,4,5,6,7])
          UHU.append(Hbond_n)
                
    elif (l_layer == 2):
      # left edge
      Hbond = get_Hbond_tensor(model,0,L,J,Fieldz,hzlist) 
      U_up  = get_umpo_up(Ulist,L,0)
      U_dn  = copy.copy(U_up)
      U_dn.permute([5,6,7,2,3,4],3)
      U_dn.setLabel([8,9,10,0,1,4])
      Hbond_n = (U_dn*Hbond)*U_up
      Hbond_n.permute([8,9,10,5,6,7],3)
      Hbond_n.setLabel([0,1,2,3,4,5])
      UHU.append(Hbond_n)

      
      for bidx in xrange(1,L-2):
        Hbond = get_Hbond_tensor(model,bidx,L,J,Fieldz,hzlist)
        U_up  = copy.copy(get_umpo_up(Ulist,L,bidx))
        U_dn  = copy.copy(U_up)
        U_dn.permute([5,6,7,8,4,2,3,20],4)
        U_dn.setLabel([9,10,11,12,4,0,1,20])
        Hbond_n = (U_dn*Hbond)*U_up
        Hbond_n.permute([9,10,11,12,5,6,7,8],4)
        Hbond_n.setLabel([0,1,2,3,4,5,6,7])
        UHU.append(Hbond_n)

      Hbond = get_Hbond_tensor(model,L-2,L,J,Fieldz,hzlist)
      U_up  = get_umpo_up(Ulist,L,L-2)
      U_dn  = copy.copy(U_up)
      U_dn.permute([5,6,7,4,2,3],3)
      U_dn.setLabel([8,9,10,4,0,1])
      Hbond_n = (U_dn*Hbond)*U_up
      Hbond_n.permute([8,9,10,5,6,7],3)
      Hbond_n.setLabel([0,1,2,3,4,5])
      UHU.append(Hbond_n)

    return UHU

       
######PAULI######PAULI######PAULI######PAULI######PAULI######PAULI
def trUP(U,P,lp):
    '''
    U  - unitary
    P  - Pauli string
    lp - length of Pauli string 
    '''
    #TODO Think if there is a way avoiding constructing Pauli strings to
    #     get the trace
    #U_ = copy.copy(U)
    #U_.setLabel([0,1,2,3])
    tr = (U*P).getBlock().sum()/(2.**lp)
    return tr

def pauli2body():
    '''
    Tensors for 1-body Pauli strings. 
    '''
    bond_dim = 2
    iden = matIden()
    sx   = matSx()
    sy   = matSy()
    sz   = matSz()
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)
    pauli1b = [iden,sx,sy,sz]
    pauli2b = []
    for i in xrange(4):
      for j in xrange(4):
        mat = uni10.otimes(pauli1b[i],pauli1b[j])
        P = uni10.UniTensor([bdi,bdi,bdo,bdo])
        P.putBlock(mat)
        P.setLabel([0,1,2,3])
        pauli2b.append(P)
    return pauli2b


def pauli3body():
    '''
    Tensors for 3-body Pauli strings. 
    '''
    bond_dim = 2
    iden = matIden()
    sx   = matSx()
    sy   = matSy()
    sz   = matSz()
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)

    pauli1b = [iden,sx,sy,sz]
    pauli3b = []
    for i in xrange(4):
      for j in xrange(4):
        for k in xrange(4):
          mat = uni10.otimes(uni10.otimes(pauli1b[i],pauli1b[j]),pauli1b[k])
          P   = uni10.UniTensor([bdi,bdi,bdi,bdo,bdo,bdo])
          P.putBlock(mat)
          pauli3b.append(P)

    return pauli3b

def pauli4body():
    '''
    Tensors for 4-body Pauli strings. 
    '''
    bond_dim = 2
    iden = matIden()
    sx   = matSx()
    sy   = matSy()
    sz   = matSz()
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)

    pauli1b = [iden,sx,sy,sz]
    pauli4b = []
    for i in xrange(4):
      for j in xrange(4):
        for k in xrange(4):
          for l in xrange(4):
            mat = uni10.otimes(uni10.otimes(uni10.otimes(pauli1b[i],pauli1b[j]),pauli1b[k]),pauli1b[l])
            P   = uni10.UniTensor([bdi,bdi,bdi,bdi,bdo,bdo,bdo,bdo])
            P.putBlock(mat)
            pauli4b.append(P)

    return pauli4b


def UHU2pauli(L,Ulist,model='Ising',J=1.0,Fieldz=1.0,hzlist=[],pauli2b=None,pauli3b=None,pauli4b=None,tol=1e-3):
    # TODO make very small terms zero so that no measurement is needed
    #      can be done when measuring too
    l_layer = len(Ulist[0])
    Nbond = L-1
    UHU = contract_UHU(L,Ulist,model,J,Fieldz,hzlist) 

    if (l_layer == 1):
      if (pauli2b is None):
        pauli2b = pauli2body()
      if (pauli4b is None):
        pauli4b = pauli4body()
      l2b = 16
      l4b = 256
      pcoefs = []
      for l in xrange(Nbond):
        if (l%2 == 0):   # EVEN BOND
          coef2 = np.zeros(l2b)
          for i in xrange(l2b):
            ans = trUP(UHU[l],pauli2b[i],2)
            if abs(ans) < tol:
              ans = 0.
            coef2[i] = ans
          pcoefs.append(coef2)
        else:            # ODD  BOND
          coef4 = np.zeros(l4b)
          for i in xrange(l4b):
            ans = trUP(UHU[l],pauli4b[i],4)
            if abs(ans) < tol:
              ans = 0.
            coef4[i] = ans
          pcoefs.append(coef4)


    elif (l_layer == 2):
      if (pauli3b is None):
        pauli3b = pauli3body()
      if (pauli4b is None):
        pauli4b = pauli4body()
      l3b = 64
      l4b = 256
      pcoefs = []

      coef3 = np.zeros(l3b)
      for i in xrange(l3b):
        ans = trUP(UHU[0],pauli3b[i],3)
        if (abs(ans)<tol):
          ans = 0.
        coef3[i] = ans
      pcoefs.append(coef3)
      for l in xrange(1,L-2):
        coef4 = np.zeros(l4b)
        for i in xrange(l4b):
          ans = trUP(UHU[l],pauli4b[i],4)
          if abs(ans) < tol:
            ans = 0.
          coef4[i] = ans
        pcoefs.append(coef4)

      coef3 = np.zeros(l3b)
      for i in xrange(l3b):
        ans = trUP(UHU[L-2],pauli3b[i],3)
        if abs(ans) < tol:
          ans = 0.
        coef3[i] = ans
      pcoefs.append(coef3)

    return pcoefs


##############Upsi##############Upsi##############Upsi##############Upsi
def UdagPsi(Ulist, psi, L):
    # TODO change this to trucated U to get exact ground state energy
    lu = L/2
    l_layer = len(Ulist[0])
    v = copy.copy(psi)


    # Applying the first layer
    label = []
    for i in xrange(lu):
        s = i*2
        U0 = copy.copy(Ulist[i][0])
        U0.permute([2,3,0,1],2)
        U0.setLabel([s+L,s+L+1,s,s+1])
        v = U0*v
        label.append([s,s+1])
    label = label[::-1]
    label = list(itertools.chain(*label))
    v.setLabel(label)
    v.permute(xrange(L),L)

    if (l_layer == 2): 
        # Applying the second layer
        label = []
        for i in xrange(lu-1):
            s=i*2+1
            U1 = copy.copy(Ulist[i][1])
            U1.permute([2,3,0,1],2)
            U1.setLabel([s+L,s+L+1,s,s+1])
            v = U1*v
            label.append([s,s+1])
        label = label[::-1]
        label = list(itertools.chain(*label))
        label = label + [0,L-1]
        v.setLabel(label)
        v.permute(xrange(L),L)
    
    return v

############MEASURE############MEASURE############MEASURE############MEASURE
def measure_pstr(Ulist,L,psi,model='Ising',J=1.0,Fieldz=1.0,hzlist=[],tol=1e-10,analys_coefs=False,outdir='./'):
    '''
    Return:
    (1) the expectation value of the Hamiltonian 
    and measurement variance under rotation:
    <psi U | U+HU| U+ psi>
    '''
    #XXX Note that I set the expectation value to be zero if
    #    the corresponding coefficient is zero.
    #UHU = contract_UHU(L, Ulist)    


    l_layer = len(Ulist[0])
    
    psiT = copy.copy(psi)
    psiT.permute(range(L),0)
    
    ave_value = [] # measurement of each Pauli string
    var_value = [] # variance of the Pauli String

    # analysis pcoefs
    
    #*************************1 Layer**************************#
    if (l_layer == 1):
        pauli2b = pauli2body()
        pauli4b = pauli4body()
        l2b = 16
        l4b = 256
        pcoefs  = UHU2pauli(L,Ulist,model,J,Fieldz,hzlist,pauli2b=pauli2b,pauli4b=pauli4b,tol=tol)
        for n in xrange(L-1):
            if (n%2==0):
                labelpsiT = range(4,5+n-1) + range(2)   + range(5+n-1,L+2)
                labelpsi  = range(4,5+n-1) + range(2,4) + range(5+n-1,L+2)
                psi.setLabel(labelpsi)
                psiT.setLabel(labelpsiT)
                coef2 = pcoefs[n]
                e2 = np.zeros(l2b)
                for i in xrange(l2b):
                    if (abs(coef2[i])<tol):
                        continue
                    else:
                        e2[i] = ((psiT*pauli2b[i])*psi).getBlock().sum()
                ave_value.append(e2)
            else:
                labelpsiT = range(8,8+n-1) + range(4)   + range(8+n-1,L+4)
                labelpsi  = range(8,8+n-1) + range(4,8) + range(8+n-1,L+4)
                psi.setLabel(labelpsi)
                psiT.setLabel(labelpsiT)
                coef4 = pcoefs[n]
                e4 = np.zeros(l4b)
                for i in xrange(l4b):
                    if (abs(coef4[i])<tol):
                        continue
                    else:
                        e4[i] = ((psiT*pauli4b[i])*psi).getBlock().sum()
                ave_value.append(e4)

    #*********************END 1 Layer**************************#
    

    #*************************2 Layers*************************#
    elif (l_layer == 2):

        pauli3b = pauli3body()
        pauli4b = pauli4body()
        pcoefs  = UHU2pauli(L,Ulist,model,J,Fieldz,hzlist,pauli3b=pauli3b,pauli4b=pauli4b,tol=tol)
        l3b = 64
        l4b = 256    
        # Left Edge
        e3 = np.zeros(l3b)
        coef3 = pcoefs[0]
        labelpsiT = range(3)   + range(6,L+3)
        labelpsi  = range(3,6) + range(6,L+3)
        psi.setLabel(labelpsi)
        psiT.setLabel(labelpsiT)
        for i in xrange(l3b):
            if (abs(coef3[i]) < tol):
                continue
            else:
                #p3 = copy.copy(pauli3b[i])
                e3[i] = ((psiT*pauli3b[i])*psi).getBlock().sum()
        ave_value.append(e3)

        # Interior Bonds
        for n in xrange(1, L-2):
            coef4 = pcoefs[n]
            e4 = np.zeros(l4b)
            labelpsiT = range(8,8+n-1) + range(4)   + range(8+n-1,L+4)
            labelpsi  = range(8,8+n-1) + range(4,8) + range(8+n-1,L+4)
            psi.setLabel(labelpsi)
            psiT.setLabel(labelpsiT)
            for i in xrange(l4b):
                if (abs(coef4[i]) < tol):
                    continue
                else:
                    e4[i] = ((psiT*pauli4b[i])*psi).getBlock().sum()
            ave_value.append(e4)

        # Right Edge
        e3 = np.zeros(l3b)
        coef3 = pcoefs[-1]
        labelpsiT = range(6,L+3) + range(3)
        labelpsi  = range(6,L+3) + range(3,6)
        psi.setLabel(labelpsi)
        psiT.setLabel(labelpsiT)
        for i in xrange(l3b):
            if (abs(coef3[i]) < tol):
                continue
            else:
                e3[i] = ((psiT*pauli3b[i])*psi).getBlock().sum()
        ave_value.append(e3)

        # Calculating the expectation values and variances
   
    #*********************END 2 Layers*************************#
    err_value = []
    for p in ave_value:
        varp = 1.-p**2.
        err_value.append(varp)
 
    
    # Save Files
    if analys_coefs:
        coef_analys(pcoefs,l_layer, w2file=True, datadir=outdir,label=model+'_L%d_lyr%d_pstr_'%(L,l_layer))

    coef_analys(err_value, l_layer, w2file=True, datadir=outdir,label='L%d_lyr%d_err_'%(L,l_layer))

    
    aveH = 0.
    errH = 0.
    for l in xrange(L-1):
        aveH += np.sum(ave_value[l]*pcoefs[l])
        errH += np.sum(err_value[l]*(pcoefs[l]**2))
    return aveH, errH

##########GROUND#STATE##########GROUND#STATE##########GROUND#STATE

def ising_ED(L,J=1.0,Fieldz=1.0):
    '''
    -J XX - h Z    (h = Fieldz)
    '''

    if (L > 12):
        raise ValueError("Exact diagonalization can only handle up to 12 sites!")
    if (L%2 != 0):
        raise ValueError("Only even number of sites is accepted!")

    npId = np.eye(2)
    npSx = np.array([0.,1.,1.,0.]).reshape(2,2)
    npSy = np.array([0.,-1.,1.,0.]).reshape(2,2)
    npSz = np.array([1.,0.,0.,-1.]).reshape(2,2)
    SxSx = np.kron(npSx,npSx)

    H = np.zeros((2**L,)*2)
    for i in xrange(0,L-1):
        H += (-1.*J)*np.kron(np.kron(np.eye(2**i),SxSx),np.eye(2**(L-i-2)))
        H += (-1.*Fieldz)*np.kron(np.kron(np.eye(2**i),npSz),np.eye(2**(L-i-1)))
    H += (-1.*Fieldz)*np.kron(np.eye(2**(L-1)),npSz) # Last Z-term

    ew,ev = np.linalg.eigh(H) # non-degenerate
    #gs = vec2uni10(ev[:,0],L)
    return ew, ev[:,0]

def heisenberg_ED(L,J=1.0,hzlist=[]):
    '''
    J(XX-YY+ZZ) - hZ (h=hzlist)
    '''

    if (L > 12):
        raise ValueError("Exact diagonalization can only handle up to 12 sites!")
    if (L%2 != 0):
        raise ValueError("Only even number of sites is accepted!")

    npId = np.eye(2)
    npSx = np.array([0.,1.,1.,0.]).reshape(2,2)
    npSy = np.array([0.,-1.,1.,0.]).reshape(2,2)
    npSz = np.array([1.,0.,0.,-1.]).reshape(2,2)
    SxSx = np.kron(npSx,npSx)
    SySy = np.kron(npSy,npSy)
    SzSz = np.kron(npSz,npSz)
    SS   = SxSx - SySy + SzSz

    H = np.zeros((2**L,)*2)
    for i in xrange(0,L-1):
        H += (1.*J)*np.kron(np.kron(np.eye(2**i),SS),np.eye(2**(L-i-2)))
        H += (-1.*hzlist[i])*np.kron(np.kron(np.eye(2**i),npSz),np.eye(2**(L-i-1)))
    H += (-1.*hzlist[-1])*np.kron(np.eye(2**(L-1)),npSz)

    ew, ev = np.linalg.eigh(H)
    return ew, ev[:,0]



def ising_variance(L,J=1.0,Fieldz=1.0):
    # TODO make it in the MPO form
    # Im in rush sorry
    sx = np.array([0.,1.,1.,0.]).reshape(2,2)
    sz = np.array([1.,0.,0.,-1.]).reshape(2,2)
    sxsx = np.kron(sx,sx)
    Egs, v = ising_ED(L,J,Fieldz)
    var = []
    for i in range(L-1):
        H_ = np.kron(np.kron(np.eye(2**i),sxsx), np.eye(2**(L-2-i)))
        var.append(1.-(np.dot(v.T,H_.dot(v)))**2)
    for i in range(L):
        H_ = np.kron(np.kron(np.eye(2**i),sz),np.eye(2**(L-1-i)))
        var.append(1.-(np.dot(v.T,H_.dot(v)))**2)

    var = np.asarray(var)
    return np.sum(var)

def heisenberg_variance(L,J=1.0,hzlist=[]):
    npId = np.eye(2)
    npSx = np.array([0.,1.,1.,0.]).reshape(2,2)
    npSy = np.array([0.,-1.,1.,0.]).reshape(2,2)
    npSz = np.array([1.,0.,0.,-1.]).reshape(2,2)
    SxSx = np.kron(npSx,npSx)
    SySy = np.kron(npSy,npSy)
    SzSz = np.kron(npSz,npSz)

    _, v = heisenberg_ED(L,J,hzlist)
    var = []
    for i in xrange(L-1):
        H_ = np.kron(np.kron(np.eye(2**i),SxSx), np.eye(2**(L-2-i)))
        var.append(1.-(np.dot(v.T,H_.dot(v)))**2)
        H_ = np.kron(np.kron(np.eye(2**i),SySy), np.eye(2**(L-2-i)))
        var.append(1.-(np.dot(v.T,H_.dot(v)))**2)
        H_ = np.kron(np.kron(np.eye(2**i),SzSz), np.eye(2**(L-2-i)))
        var.append(1.-(np.dot(v.T,H_.dot(v)))**2)

    for i in xrange(L):
        H_ = np.kron(np.kron(np.eye(2**i),npSz),np.eye(2**(L-1-i)))
        var.append( (hzlist[i]**2)*(1.-(np.dot(v.T,H_.dot(v)))**2) )

    var = np.asarray(var)
    return np.sum(var)
         
        
#########UTIL#########UTIL#########UTIL#########UTIL#########UTIL#########UTIL
def str2int(mystr, base):
    l = len(mystr)
    str_n = mystr[::-1]
    myint = 0
    for i in xrange(l):
        myint += str_n[i]*(base**i)
    return int(myint)

def coef_analys(pcoefs,l_layer,w2file=False, datadir='./', label=''):

    lc = len(pcoefs)

    if (l_layer == 1):
        zlabel = [0,3]
        nzlabel = [1,2]
        zcoef = []
        nzcoef = []
        for p in xrange(lc):
            zcoef.append(pcoefs[p][zlabel])
            nzcoef.append(pcoefs[p][nzlabel])

    elif (l_layer == 2):
        l3b = 64
        l4b = 256    
        zlabel3 = []
        zlabel4 = []
        for i in [0,3]:
          for j in [0,3]:
            for k in [0,3]:
              zlabel3.append(str2int([i,j,k],4))
              for l in [0,3]:
                zlabel4.append(str2int([i,j,k,l],4))
        nzlabel3 = range(l3b)
        nzlabel4 = range(l4b)
        for i in zlabel3:
            nzlabel3.remove(i)
        for i in zlabel4:
            nzlabel4.remove(i)
        zcoef  = []
        nzcoef = []

        zcoef.append(pcoefs[0][zlabel3])
        nzcoef.append(pcoefs[0][nzlabel3])

        for p in xrange(1,lc-1):
            zcoef.append(pcoefs[p][zlabel4])
            nzcoef.append(pcoefs[p][nzlabel4])
        zcoef.append(pcoefs[-1][zlabel3])
        nzcoef.append(pcoefs[-1][nzlabel3])
    
    if w2file:
        with open (datadir+'/'+label+'_Z.txt', 'wb') as zfp:
            pickle.dump(zcoef,zfp)
        with open (datadir+'/'+label+'_NZ.txt', 'wb') as nzfp:
            pickle.dump(nzcoef,nzfp)
    else:
        return zcoef, nzcoef 

def uni2mat(uni,lmat):
    # uni10 matrix to numpy array
    try:
        M = uni.getBlock()
    except:
        M = uni
    npmat = np.zeros((lmat,lmat))
    for i in xrange(lmat):
        for j in xrange(lmat):
            npmat[i,j] = M[i*lmat+j]

    return npmat

def vec2uni10(v,L):
    bond_dim = 2
    vl = list(v)
    vl10 = uni10.Matrix(2**L,1,vl)
    bdi = uni10.Bond(uni10.BD_IN,  bond_dim)
    bdo = uni10.Bond(uni10.BD_OUT, bond_dim)
    vtensor = uni10.UniTensor([bdi]*L,'GS')
    vtensor.putBlock(vl10)
    return vtensor

def get_uni10diag(uni,lmat):
    try:
        M = uni.getBlock()
    except:
        M = uni
    npdiag = np.zeros(lmat)
    for i in xrange(lmat):
        npdiag[i] = M[i*(lmat+1)]
    return npdiag
        

#########SPECTRUM#########SPECTRUM#########SPECTRUM#########SPECTRUM
def get_spectrum(Ulist,model,L,UHU=None,J=1.,Fieldz=1.,hzlist=[],sort_spec=True):
    # Get the diagonal of rotated Hamiltonian U+HU
    if UHU is None:
        UHU = contract_UHU(L, Ulist, model, J, Fieldz, hzlist)
    diag = np.zeros(2**L) 
    d0 = get_uni10diag(UHU[0],8)
    diag += np.kron(d0,np.ones(2**(L-3)))
    for b in xrange(1,L-2):
        d0 = get_uni10diag(UHU[b],16)
        diag += np.kron(np.ones(2**(b-1)), np.kron(d0, np.ones(2**(L-b-3))))
    d0 = get_uni10diag(UHU[-1],8)
    diag += np.kron(np.ones(2**(L-3)),d0)
    if sort_spec:
        diag = np.sort(diag)
    
    return diag
    
