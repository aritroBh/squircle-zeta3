import os, pickle
import sympy as sp, mpmath as mp
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/coef_acc_raw.pkl','rb') as f: coef_acc = pickle.load(f)
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
# take one coefficient
mk = list(coef_acc.keys())[0]
co = coef_acc[mk]
print("monomial:", mk)
print("coef free syms:", sorted([str(s) for s in co.free_symbols]))
print("lA0 in free syms?", lA0 in co.free_symbols)
print("type of one lA0 from free_symbols:", [type(s) for s in co.free_symbols if str(s)=='lA0'])
# test subs
test = co.subs({lA0: sp.Float(-0.5), lA1: sp.Float(-1.0), lB: sp.Float(-1.0)})
print("after log subs, free syms:", sorted([str(s) for s in test.free_symbols]))
