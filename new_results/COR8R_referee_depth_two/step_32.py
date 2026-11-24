import cor8r_sym as S
import sympy as sp
from fractions import Fraction
import itertools
cc = S.cc
A = 1/(1+cc); B = 2/(1+cc)

# reproduce the failing system for length 2 over alphabet {0,A,B}
L = 2
values = [A,B]
unk = S._bad_words(L, values)
print("unknowns:", unk)
alphabet = [sp.Integer(0)]+values
for lu in range(1,L):
    lv = L-lu
    for u in itertools.product(alphabet, repeat=lu):
        su = S.symG(u)
        for v in itertools.product(alphabet, repeat=lv):
            sv = S.symG(v)
            prod = S.S_shuffle(su,sv)
            print("u=",u,"v=",v," interleavings:", S._shuffle_words(u,v))
