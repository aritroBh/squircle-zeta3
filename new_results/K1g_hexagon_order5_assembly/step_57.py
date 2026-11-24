import os, pickle, time
import sympy as sp
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/coef_acc_raw.pkl','rb') as f: coef_acc = pickle.load(f)
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
# Constraint: A0 = 1 - B - B*kap
constr = {A0: 1 - B - B*kap}
t0=time.time()
print("Checking pure-rational & pure-z pole coefficients vanish under A0=1-B-B*kap:")
nchecked=0; nzero=0
results={}
for mk in sorted(coef_acc.keys(), key=lambda m:m[1]):
    at,xe,ze = mk
    if at==() and xe<0:  # pure (rational or z) poles
        co = coef_acc[mk]
        co_sub = sp.simplify(co.subs(constr))
        is_zero = (co_sub==0)
        results[mk]=is_zero
        nchecked+=1
        if is_zero: nzero+=1
        print(f"  x^({xe}) z{ze}: {'ZERO' if is_zero else 'NONZERO -> '+str(sp.factor(co_sub))[:80]}")
print(f"\n{nchecked} pure pole coefficients checked, {nzero} vanish. time {time.time()-t0:.0f}s")
