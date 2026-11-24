symsrc = '''"""Weight-2 and weight-3 symbols of Acc expressions, via q1s_symbol.sym_I on words
and shuffle with lv-monomial symbols. Splitting coefficients as A + B sqrt2 (+C i + D i sqrt2).
"""
from fractions import Fraction as Fr
from q1s_exact import *
from q1s_symbol import sym_I, letter_vec, tadd, tscale, t_tensor, t_shuffle, t_conj, teq, tsub

def split_sqrt2(c):
    """F-element c = A + B sqrt2 + C i + D i sqrt2 -> (A,B,C,D) Fractions."""
    x0,x1,x2,x3=c
    return (x0, (x1-x3)/2, x2, (x1+x3)/2)

def lv_sym(f):
    """symbol (vector) of ('lv',x) = log(x)."""
    return letter_vec(f[1])

def monomial_sym(lvm):
    T={():Fr(1)}
    for f in lvm:
        T=t_shuffle(T, {(k,):v for k,v in lv_sym(f).items()})
    return T

def acc_symbol(acc, weight):
    """returns dict (A,B,C,D) of tensors (dict tuple->Fr) at given weight (2 or 3)."""
    out=[{},{},{},{}]
    def add(T, c):
        A,B,C,D=split_sqrt2(c)
        for i,part in enumerate((A,B,C,D)):
            if part:
                out[i]=tadd(out[i], tscale(T, part))
    for (w,lvm),c in acc.W.items():
        if len(w)+len(lvm)!=weight: continue
        Sw=sym_I(w, ZERO, ONE)
        Sm=monomial_sym(lvm)
        add(t_shuffle(Sw,Sm), c)
    for lvm,c in acc.P.items():
        if len(lvm)!=weight: continue
        add(monomial_sym(lvm), c)
    return out
'''
open('q3_symbol.py','w').write(symsrc)
import importlib, q3_engine, q1s_exact, q3_stage2, q3_symbol, q1s_symbol
for m in (q1s_exact, q1s_symbol, q3_engine, q3_stage2, q3_symbol): importlib.reload(m)
from q3_stage2 import build_V1, build_V3p
from q3_symbol import acc_symbol
from q1s_symbol import t_conj, teq
import pickle
V1=build_V1(); V3=build_V3p()
t0=__import__('time').time()
S3_V1 = acc_symbol(V1, 3)
S3_V3 = acc_symbol(V3, 3)
print("time:", __import__('time').time()-t0)
for nm,S in [("V1",S3_V1),("V3'",S3_V3)]:
    A,B,C,D = S
    print(nm, ": |A|=%d |B|=%d |C|=%d |D|=%d  conj(A)==A:%s conj(B)==B:%s"%(
        len(A),len(B),len(C),len(D), teq(t_conj(A),A), teq(t_conj(B),B)))
pickle.dump({'S3_V1':S3_V1,'S3_V3':S3_V3}, open('q3_symbols.pkl','wb'))
