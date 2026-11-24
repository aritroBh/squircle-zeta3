import os, pickle
import sympy as sp
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f, 'E5_polys.pkl'), 'rb') as f:
    E5_polys = pickle.load(f)
with open(os.path.join(k1f, 'jtab5.pkl'), 'rb') as f:
    jtab5 = pickle.load(f)
with open(os.path.join(k1f, 'jloc_sym.pkl'), 'rb') as f:
    jloc = pickle.load(f)
E5 = E5_polys[5]

# --- Symbols ---
# Edge data
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
# atom argument
x = sp.Symbol('x')
# HPL atoms (13 Lyndon) - use distinct names
atom_names = ['H1','Li2','Li3','H011','Li4','H0011','H0111','Li5',
              'H00011','H00101','H00111','H01011','H01111']
atoms = {nm: sp.Symbol(nm) for nm in atom_names}
z2,z3,z5 = sp.symbols('z2 z3 z5')

# namespace for parsing jtab5 (non-local): atoms + x + z2,z3,z5
ns_nl = dict(atoms)
ns_nl['x'] = x
ns_nl['z2']=z2; ns_nl['z3']=z3; ns_nl['z5']=z5

# parse jtab5
J5 = {}
for k, s in jtab5.items():
    J5[k] = sp.sympify(s, locals=ns_nl)
print("parsed jtab5, sample:")
print("J5[(0,0,0,1)] =", J5[(0,0,0,1)])
print("J5[(3,1,1,1)] free syms:", sorted([str(s) for s in J5[(3,1,1,1)].free_symbols]))

# namespace for jloc (local): DEN, HH1..HH5, z2,z3,z5
DEN = sp.Symbol('DEN')
HH = {i: sp.Symbol(f'HH{i}') for i in range(1,6)}
ns_loc = {'DEN':DEN, 'z2':z2,'z3':z3,'z5':z5}
for i in range(1,6):
    ns_loc[f'H{i}'] = HH[i]
Jloc = {}
for k, s in jloc.items():
    Jloc[k] = sp.sympify(s, locals=ns_loc)
print("\nparsed jloc, sample:")
print("Jloc[(1,1)] =", Jloc[(1,1)])
print("free syms:", sorted([str(s) for s in Jloc[(1,1)].free_symbols]))
