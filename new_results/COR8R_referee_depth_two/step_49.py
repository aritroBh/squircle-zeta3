src = open('cor8r_sym.py').read()
old = """    num, den = sp.fraction(f)
    vec = [0]*NB"""
new = """    num, den = sp.fraction(f)
    vec = [0]*5   # slot 0 is [-1] (killed on return)"""
assert old in src; src = src.replace(old, new)
open('cor8r_sym.py','w').write(src)
# quick test
import importlib, cor8r_sym as S
importlib.reload(S)
import sympy as sp
cc = S.cc
print(S.factor_Qc(2))
print(S.factor_Qc(-cc))          # [-c] = [c] since [-1]=0
print(S.factor_Qc((1-cc)/(2*cc)))
print(S.factor_Qc(2/(1+cc)))
print(S.factor_Qc((cc-1)/(1+cc)))
