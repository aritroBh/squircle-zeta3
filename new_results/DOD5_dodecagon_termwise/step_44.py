import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from fractions import Fraction as Fr
import sympy as sp
from sympy import Rational, sqrt, symbols, Function, Add
from dod3_R import *   # build_R, TRIPLES, logfac etc (lu = log(2-sqrt3) convention)
t0=time.time()
base2 = "<project>/research/08_theory/"
ac = pickle.load(open(base2+"DOD2_dodecagon_finish/atom_coefs_pq.pkl","rb"))

# rebuild sum c^T3a * R_sym exactly (DOD3 step_38 logic)
li3_content = {}
li2log_part = 0; logcube_part = 0
for (i,j,k) in TRIPLES:
    p,q = ac[('T3a',(i,j,k))]
    c = Rational(p) + Rational(q)*r3
    R = build_R(i,j,k)
    cR = sp.expand(c*R)
    for term in Add.make_args(cR):
        a3 = term.atoms(Li3A); a2 = term.atoms(Li2A)
        if a3:
            assert len(a3)==1
            pt = list(a3)[0].args[0]
            coef = sp.simplify(term / list(a3)[0])
            ppt = pair_of(pt); cpq = pair_of(coef)
            li3_content[ppt] = padd(li3_content.get(ppt, P(0,0)), cpq)
        elif a2:
            li2log_part += term
        else:
            logcube_part += term
li2log_part = sp.expand(li2log_part); logcube_part = sp.expand(logcube_part)
print("rebuilt sum cR: t=%.0fs"%(time.time()-t0))
# compare with DOD3's pickle
R3c = pickle.load(open(base2+"DOD3_dodecagon_symbol/R_li3_content.pkl","rb"))
mine = {k:v for k,v in li3_content.items() if not peq(v,P(0,0))}
theirs = {k:v for k,v in R3c['li3_content'].items() if not peq(v,P(0,0))}
print("li3 points: mine", len(mine), " theirs", len(theirs), " identical:", all(peq(mine[k],theirs.get(k,P(9,9))) for k in mine) and set(mine)==set(theirs))
print("li2log identical:", sp.simplify(li2log_part - sp.sympify(R3c['li2log']))==0)
print("logcube identical:", sp.simplify(logcube_part - sp.sympify(R3c['logcube']))==0)
print("\nR-piece Li3 content (12 points max):")
for pt,c in sorted(mine.items(), key=lambda kv: pfloat(kv[0])):
    print("   Li3(%-14s ~%+.8f)  coeff = (%s) + (%s)*sqrt3"%(str(pt), pfloat(pt), psplit(c)[0], psplit(c)[1]))
pickle.dump({'li3_content':mine,'li2log':str(li2log_part),'logcube':str(logcube_part)}, open('dod5_R_content.pkl','wb'))
print("saved dod5_R_content.pkl, elapsed %.0fs"%(time.time()-t0))
