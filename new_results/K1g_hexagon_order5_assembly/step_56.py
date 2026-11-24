import os, pickle
import sympy as sp
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/coef_acc_raw.pkl','rb') as f: coef_acc = pickle.load(f)
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
# Look at pure-rational monomials (no atoms, no z) with xe<0 -- these are genuine poles
print("Pure monomials (atoms=(), ze=(0,0,0)) with x-poles:")
for mk in sorted(coef_acc.keys(), key=lambda m:(m[1])):
    at,xe,ze = mk
    if at==() and ze==(0,0,0) and xe<0:
        co = coef_acc[mk]
        print(f"  x^({xe}): coef has {len(sp.Add.make_args(co))} terms")
# Also pure z monomials with poles
print("\nPure z-monomials (atoms=(), ze!=(0,0,0)) with x-poles:")
found=False
for mk in sorted(coef_acc.keys()):
    at,xe,ze = mk
    if at==() and ze!=(0,0,0) and xe<0:
        print(f"  x^({xe}) z{ze}: present")
        found=True
if not found: print("  (none)")
# Count monomials by (has_atoms, xe sign)
import collections
cc=collections.Counter()
for mk in coef_acc:
    at,xe,ze=mk
    cc[(len(at)>0, xe<0)] += 1
print("\nmonomial count by (has_atoms, has_pole):", dict(cc))
