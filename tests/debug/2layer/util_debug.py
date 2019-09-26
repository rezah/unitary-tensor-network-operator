import numpy as np
import sys
sys.path.append("/Users/sunchong/CodingDir/TNMBL")
from qsampl import *
from pauli2body import *
import pyUni10 as uni10
from scipy import linalg as sla
import copy 


def test_getuni10diag():
    M = uni10.Matrix(3,3,[1.,2.,3.,4.,5,6,7,8,9])
    x = get_uni10diag(M,3) 
    print x

def test_p1body():
    p1b = pauli1body()
    for p in p1b:
        print p

test_p1body()
exit()


test_getuni10diag()
    



