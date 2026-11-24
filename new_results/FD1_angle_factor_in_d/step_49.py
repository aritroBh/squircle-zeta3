import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
from fractions import Fraction as Fr
from math import factorial
data=pickle.load(open(fd1+'/kern_sq_11_20.pkl','rb'))['kernels']
K=data
# inspect format
sample = K[11][(10,1)] if (10,1) in K[11] else K[11][(11,0)]
print(type(sample), sample if isinstance(sample,dict) else sample.keys())
key0=list(K[11].keys())[0]
P=K[11][key0]
print("key", key0, "parts:", list(P.keys()) if isinstance(P,dict) else type(P))
print("real sample:", list(P['real'].items())[:2])
