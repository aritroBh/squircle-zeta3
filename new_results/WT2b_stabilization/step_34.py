import sys, pickle
sys.path.insert(0,'.')
from fastred import *
record = pickle.load(open('record_W2b.pkl','rb'))
from fractions import Fraction as Fr

# top-order +i coefficient of HigherPart of the INITIAL data R_ij = W_ij at k=i+j
def rho(p,q):
    d=p+q
    rec = record[(d,(p,q))]
    H = rec['H'][1]
    top = 2*d-2
    return H[top-1] if len(H)>=top else ZERO

# my derived formula: rho_pq = K_d (-1)^p [C(d,p)+C(d,q)] ; fit K_d from data
vals={}
for d in range(2,8):
    for p in range(d,(d-1)//2,-1):
        q=d-p
        vals[(p,q)] = rho(p,q)
        print(f"rho_{p}{q} = {vals[(p,q)].sympy()}")
print()
# ratio test: K_d from (d,0):  rho_d0 = K_d (-1)^d [1+1] = 2 K_d (-1)^d
for d in range(2,8):
    r = vals[(d,0)]
    Kd = r / (2 * ((-ONE)**d))
    print(f"d={d}: K_d = {Kd.sympy()}    -i(2d-1)/(d+1) ratio vs next:")
# check ratio K_{d+1}/K_d vs -i(2d-1)/(d+1) and vs -i(2d-1)/(2(d+1))
print()
IQd = QSI(0,0,1,0)
Klist = {}
for d in range(2,8):
    Klist[d] = vals[(d,0)] / (2*((-ONE)**d))
for d in range(2,7):
    rat = Klist[d+1]/Klist[d]
    pred1 = QSI(0,0,-(2*d-1),0)/QSI(d+1)
    pred2 = QSI(0,0,-(2*d-1),0)/QSI(2*(d+1))
    print(f"d={d}: K_{{d+1}}/K_d = {rat.sympy()}   == -i(2d-1)/(d+1)? {rat==pred1}   == -i(2d-1)/(2(d+1))? {rat==pred2}")
