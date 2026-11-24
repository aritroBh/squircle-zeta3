import pickle, sympy as sp, sys
d,m=sp.symbols('d m')
loc={'d':d,'m':m,'I':sp.I}
A=pickle.load(open('lib/regen_F78.pkl','rb'))
B=pickle.load(open('orig_cache/forced_F78.orig.pkl','rb'))
print('same keys', set(A)==set(B), len(A))
print('byte-identical strings:', all(A[k]==B[k] for k in A))
bad=[]
for k in sorted(A):
    a=sp.sympify(A[k],locals=loc); b=sp.sympify(B[k],locals=loc)
    z=sp.cancel(sp.together(a-b))
    if z!=0: bad.append(k)
print('symbolic difference nonzero at:', bad)
