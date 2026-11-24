import pickle, sympy as sp, sys, os
wd = '<project>/research/08_theory/TOW5R_review'
sys.path.insert(0, wd+'/tow5_copy')
a9 = pickle.load(open(wd+'/tow5_copy/a9_closed.pkl','rb'))
a10 = pickle.load(open(wd+'/tow5_copy/a10_closed.pkl','rb'))
print('a9_closed.pkl keys:', list(a9.keys()) if isinstance(a9, dict) else type(a9))
print('a10_closed.pkl keys:', list(a10.keys()) if isinstance(a10, dict) else type(a10))
d, m = sp.symbols('d m')
import tow5lib as T5
A9 = T5.arho(9); A10 = T5.arho(10)
for k, v in (a9.items() if isinstance(a9, dict) else []):
    print('a9 item:', k, str(v)[:120])
for k, v in (a10.items() if isinstance(a10, dict) else []):
    print('a10 item:', k, str(v)[:120])
