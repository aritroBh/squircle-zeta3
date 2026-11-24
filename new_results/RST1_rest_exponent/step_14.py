"""RST1 step_14: validate v2 assembly; monotonicity of every normalized piece on grid."""
import mpmath as mp
from rst1_common2 import *
rows=load_mk()
print(" k      C(k)=total bound/target     |rest|/(C*target)")
for k in [99,150,250,400]:
    C=totalC(k); mm=Mmain(k); rest=abs(rows[k]-mm)
    print(f"{k:4d}   {mp.nstr(C,10):>12s}        {mp.nstr(rest/(C*tg(k)),8)}")
print()
print("piece values at k=99:", {n: mp.nstr(f(99),8) for n,f in [('q2',P_q2),('mid',P_mid),('mid2',P_mid2),('tail',P_tail),('T2',P_T2),('T3J',P_T3J)]})
print("total C at 99:", mp.nstr(totalC(99),10))
# monotonicity grid check for each piece
import numpy as np
for nm,f in [('q2',P_q2),('mid',P_mid),('T2',P_T2)]:
    bad=0
    for k in range(99,1200):
        if f(k+1)>f(k)*(1+1e-11): bad+=1
    print(f"monotonicity violations {nm}: {bad}")
for nm,f in [('mid2',P_mid2),('tail',P_tail),('T3J',P_T3J)]:
    bad=0
    for k in range(99,800,2):
        if f(k+2)>f(k)*(1+1e-11): bad+=1
    print(f"monotonicity violations {nm}: {bad}")
print("large-k C:", [(K, mp.nstr(totalC(K),8)) for K in [99,400,2000,20000]])
