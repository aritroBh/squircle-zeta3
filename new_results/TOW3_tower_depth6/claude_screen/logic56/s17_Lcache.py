import sys, pickle
sys.dont_write_bytecode = True
import sympy as sp
import tow3lib as T3
d, m = T3.d, T3.m
Ls = pickle.load(open('L_simplified.pkl', 'rb'))
Lp = pickle.load(open('L_parsed.pkl', 'rb'))
print('L_simplified.pkl:', type(Ls), sorted(Ls.keys())[:10], type(Ls[1]))
print('L_parsed.pkl    :', type(Lp), sorted(Lp.keys())[:10], type(Lp[1]))
loc = {'d': d, 'm': m, 'I': sp.I}
ok = {}
for r in range(1, 9):
    a = Ls[r]; a = sp.sympify(a, locals=loc) if isinstance(a, str) else a
    b = Lp[r]; b = sp.sympify(b, locals=loc) if isinstance(b, str) else b
    ok[r] = sp.cancel(sp.together(a-b)) == 0
print('L_simplified == L_parsed for r=1..8:', ok)
