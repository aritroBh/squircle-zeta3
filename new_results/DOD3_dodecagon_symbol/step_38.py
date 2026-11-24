import sys, os, pickle, time
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from dod3_R import *

base = os.path.dirname(os.getcwd())
d2p = os.path.join(base,'DOD2_dodecagon_finish')
ac = pickle.load(open(os.path.join(d2p,'atom_coefs_pq.pkl'),'rb'))

t0=time.time()
# build sum over 12 triples of cT3a * R_sym, collect Li3A atoms with exact F-coefficients
li3_content = {}   # pair point -> F-coeff pair (m,n,d) meaning (m+n r3)/d
li2log_part = 0    # sympy expr (products): Li2A*log terms
logcube_part = 0   # sympy expr: pure log cubes
for (i,j,k) in TRIPLES:
    p,q = ac[('T3a',(i,j,k))]
    c = Rational(p) + Rational(q)*r3
    R = build_R(i,j,k)
    cR = sp.expand(c*R)
    # split cR into Li3A terms, Li2A*log terms, logcube terms
    for term in Add.make_args(cR):
        a3 = term.atoms(Li3A)
        a2 = term.atoms(Li2A)
        if a3:
            assert len(a3)==1
            pt = list(a3)[0].args[0]
            coef = sp.simplify(term / list(a3)[0])
            ppt = pair_of(pt)
            cpq = pair_of(coef)
            if ppt in li3_content:
                li3_content[ppt] = padd(li3_content[ppt], cpq)
            else:
                li3_content[ppt] = cpq
        elif a2:
            li2log_part += term
        else:
            logcube_part += term
print("built sum cT3a*R_sym: t=%.0fs"%(time.time()-t0))
print("\nR-piece Li3 content (points + exact F-coefficients):")
tot_terms = 0
for pt,c in sorted(li3_content.items(), key=lambda kv: pfloat(kv[0])):
    if peq(c,P(0,0)): continue
    tot_terms += 1
    p_,q_ = psplit(c)
    print("   Li3(%-14s ~%+.8f)  coeff = %s + %s*sqrt3   [pair %s]"%(str(pt), pfloat(pt), p_, q_, c))
print("distinct Li3 points in R (nonzero coeff):", tot_terms)
# also check zero-coeff points
print("\npoints with zero total coeff:", [str(pt) for pt,c in li3_content.items() if peq(c,P(0,0))])
print("\nLi2*log part #terms:", len(Add.make_args(sp.expand(li2log_part))))
print("logcube part #terms:", len(Add.make_args(sp.expand(logcube_part))))
pickle.dump({'li3_content':li3_content,
             'li2log':str(li2log_part),'logcube':str(logcube_part)}, open('R_li3_content.pkl','wb'))
print("saved R_li3_content.pkl")
