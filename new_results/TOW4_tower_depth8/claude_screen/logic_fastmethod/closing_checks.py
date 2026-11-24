import sys, random, pickle
import sympy as sp
import pw
from pw import G, Fu, TAU1, TAU2, rand_pt, to_G
from fractions import Fraction as Fr
rng = random.Random(777)
d,m = pw.d_, pw.m_
loc={'d':d,'m':m,'I':sp.I}
cache = pickle.load(open('lib/forced_F78.pkl','rb'))
# 1. pointwise engine == cached symbolic F at random points
bad=0; n=0
for _ in range(4):
    D,M = rand_pt(rng)
    Dsp = sp.Rational(D.a.numerator, D.a.denominator); Msp = sp.Rational(M.a.numerator, M.a.denominator)
    for (delta,s),txt in cache.items():
        ex = sp.sympify(txt, locals=loc)
        v = to_G(ex.subs({d:Dsp, m:Msp}))
        w = Fu(delta,s,D,M)
        n+=1
        if not (v-w).iszero(): bad+=1; print('MISMATCH',delta,s)
print(f'cached symbolic F vs independent pointwise F: {n} (delta,s,point) cells, mismatches {bad}')

# 2. mutation tests of the closing conditions: perturb one input, closing must fail
def closing_vals(npt=3):
    out=[]
    for _ in range(npt):
        pw._F.clear(); pw._acache.clear()
        D,M = rand_pt(rng)
        out.append((Fu(7,4,D,M).iszero(), Fu(8,5,D,M).iszero(), Fu(6,4,D,M).iszero(), Fu(5,3,D,M).iszero()))
    return out
print('baseline closing zero? (F74,F85,F64,F53):', closing_vals())
import copy
t1=copy.deepcopy(pw.TAU1); t2=copy.deepcopy(pw.TAU2)
for name,delta_ in [('tau1[8] += 1/128',('t1',8)),('tau2[8] += 1/128',('t2',8)),('tau1[7] += 1/64',('t1',7)),('tau2[6] += 1/64',('t2',6)),('tau1[3] += 1/64',('t1',3))]:
    pw.TAU1.update(copy.deepcopy(t1)); pw.TAU2.update(copy.deepcopy(t2))
    which,n_ = delta_
    tab = pw.TAU1 if which=='t1' else pw.TAU2
    tab[n_] = tab[n_] + G(Fr(1,128))
    print('MUTANT', name, ' closing zero? (F74,F85,F64,F53):', closing_vals())
pw.TAU1.update(t1); pw.TAU2.update(t2)
# mutate L cache entries
orig = dict(pw._Lpoly)
for r in (8,7,6,5):
    pw._Lpoly.clear(); pw._Lpoly.update(orig)
    num,den = orig[r]
    mon0 = ((0,0), G(1))   # add 1*(den) to numerator => L_r += 1
    pw._Lpoly[r] = (num + [((mon0[0]), mon0[1]*G(0))] + [ (mn, c) for mn,c in den ], den)  # numerator + den  => L_r + 1
    print('MUTANT L_%d += 1'%r, ' closing zero? (F74,F85,F64,F53):', closing_vals())
pw._Lpoly.clear(); pw._Lpoly.update(orig)
print('restored baseline:', closing_vals(2))
