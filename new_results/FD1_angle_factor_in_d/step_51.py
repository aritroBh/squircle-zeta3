import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
from fractions import Fraction as Fr
K=pickle.load(open(fd1+'/kern_sq_11_20.pkl','rb'))['kernels']
poles=set()
for k in K:
    for key in K[k]:
        for rho in K[k][key]['real']:
            poles.add(rho)
from mpmath import sqrt
def val(rho):
    return float(Fr(*rho[0])) + float(Fr(*rho[1]))*float(sqrt(2))
for p in sorted(poles, key=val):
    print(p, " ~ ", val(p))
