import numpy as np
import sys
sys.path.append('/Users/sunchong/CodingDir/TNMBL/')
from utils import *
import copy

#mystr = sys.argv[1]
#print str2pauli(mystr)

#coef = np.random.rand(16)
#print print_paulistr_expansion(coef,tol=1e-5)

#mat = pauli_multiply_dict()
#print mat[0]
#print mat[1]
#
#p1 = sys.argv[1]
#p2 = sys.argv[2]
#s,n = multiply_pstr(p1,p2)
#print s
#print n
L = 6
v = np.random.rand(2**L)
v = v/np.linalg.norm(v)
maxm = 10
mps = vec2mps(v,L,maxm)
print mps
#print mps[0].getBlock()


