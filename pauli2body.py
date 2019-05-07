import pyUni10 as uni10
import copy
import numpy as np

'''
Functions to translate (2-body) unitary blocks into Pauli strings and 
do the contraction with Pauli strings.
'''

######PAULI######PAULI######PAULI######PAULI######PAULI######PAULI
def pauli2body():
    # TODO figure out an easier way to print the labels
    '''
    Tensors for 2-body Pauli strings. 
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
    pauli2b_labels = [['II','IX','IY','IZ'],\
                      ['XI','XX','XY','XZ'],\
                      ['YI','YX','YY','YZ'],\
                      ['ZI','ZX','ZY','ZZ']]            
    for i in range (4):
        for j in range(4):
            mat = uni10.otimes(pauli1b[i],pauli1b[j])
            P   = uni10.UniTensor([bdi,bdi,bdo,bdo],pauli2b_labels[i][j])
            P.putBlock(mat)
            pauli2b.append(P) 

    return pauli2b, pauli2b_labels

def pauli_table():
    '''
    Construct the Pauli multiplication table.
    Return a 4x4 signed Matrix. 
    +/- index.
    '''
    table = np.zeros((4,4),dtype=np.int32)
    for i in range(4):
        table[i,i] = 0
        table[0,i] = i
    table[1,2] = 3
    table[1,3] = 2
    table[2,3] = 1
    table = table + table.T

    table[2,1] *= -1
    table[3,1] *= -1
    table[3,2] *= -1
    
    return table
            
    
def unitary2pauli(U,P,lp):
    '''
    U  - unitary
    P  - Pauli string
    lp - length of Pauli string 
    '''
    U_ = copy.copy(U)
    U_.setLabel([0,1,2,3])
    trUP = (U_*P).getBlock().sum()/(2.*lp)
    return trUP

def ulist2pauli(Ulist,lp):
    '''
    Turn the Ulist into a list of vectors.
    '''
    #TODO Call pauli2body() in a higher fuction to reduce the number of calls
    #TODO Expecting most coefficients zero, so we can use dictionary to take
    #     down only non-zero coefficients and their indices.
    ptensor, plabel = pauli2body() 
    l0 = len(Ulist)
    l1 = len(Ulist[0])
    np2 = len(ptensor)
    print np2
    
    Ucoef = []
    
    if l1 != 2:
        raise Exception("Only 2-layer U_list is accepted!") 
    for i in range(l0):
        ucoef_ = []
        for j in range(l1):
            coef = np.zeros(np2)
            U = Ulist[i][j]
            for k in range(np2):
                coef[k] = unitary2pauli(U,ptensor[k],lp)
                if (abs(coef[k])<1e-10):
                    coef[k] = 0.   # setting small coefficients to zero
            ucoef_.append(coef)
        Ucoef.append(ucoef_)     

    return Ucoef

def pauli_str_contract(strA,strB,posA,posB,table=None):
    '''
    Only continuous contractions are allowed.
    A is in front of B.
    '''
    if table is None:
        table = pauli_table() 
    if len(posA) != len(posB):
        raise Exception("Contraction indices of A and B should have same size!")
    lc = len(posA)
    ls = len(posA)+len(posB)-len(posA)
    newstr = ''
    coef   = 1.
    if (posA[0] > posB[0]):
        newstr += strA[:posA[0]]
        for i in range(lc):
            pA = int(strA[posA[i]])
            pB = int(strB[posB[i]])
            newstr += str(abs(table[pA,pB]))
            sgn     = np.sign(table[pA,pB])
            if abs(sgn) < 1e-10:
                sgn = 1.
            coef   *= sgn
        newstr += strB[(posB[-1]+1):]

    else:
        newstr += strB[:posB[0]]
        for i in range(lc):
            pA = int(strA[posA[i]])
            pB = int(strB[posB[i]])
            newstr += str(abs(table[pA,pB]))
            sgn     = np.sign(table[pA,pB])
            if abs(sgn) < 1e-10:
                sgn = 1.
            coef   *= sgn
        newstr += strA[(posA[-1]+1):]

    return newstr, coef
        
    
 
