import os, pickle, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
mp.mp.dps=50
import sympy as sp
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
res=pickle.load(open('gamma_trace2.pkl','rb'))
ks=sorted(res.keys())
print("keys sample:", ks[:3], "...", [k for k in ks if 0.008<k<0.02], [k for k in ks if 1.5<k<2], [k for k in ks if 250<k<350])
