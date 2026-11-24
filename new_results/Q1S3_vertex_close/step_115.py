import pickle, itertools
from fractions import Fraction as Fr
import sympy as sp
from q1s_symbol import tadd, tscale, t_conj, t_shuffle, t_from_vec, t_tensor, letter_vec
from q1s_exact import BASISEL, fconj, feq
B4=('a','p2','q1','q2')
sols=pickle.load(open('sunit_li3.pkl','rb'))
sol=pickle.load(open('q3_solve_min.pkl','rb'))['V3p']
S=pickle.load(open('q3_symbols.pkl','rb'))
# rebuild generator symbols
def li3_sym(pr):
    i,j=pr
    return tadd(sols[i]['S'], t_conj(sols[i]['S']))
def dec_sym(d):
    u,v,w=d
    t={}
    for perm in [(u,v,w),(v,u,w),(v,w,u)]: t[perm]=t.get(perm,Fr(0))+1
    return t
tot={}
for (k,v),c in sol:
    c=Fr(c)
    g = li3_sym(v) if k=='li3' else dec_sym(v)
    tot=tadd(tot, tscale(g,c))
tgt=S['S3_V3'][0]
keys=set(tot)|set(tgt)
diff={k:tot.get(k,Fr(0))-tgt.get(k,Fr(0)) for k in keys}
print("symbol identity exact:", all(c==0 for c in diff.values()), " max dev:", max(abs(c) for c in diff.values()))
