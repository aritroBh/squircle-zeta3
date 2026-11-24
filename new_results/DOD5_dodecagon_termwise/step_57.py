import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
t0=time.time()
mp.mp.dps = 620
red = pickle.load(open('dod5_all_li3_reductions.pkl','rb'))
xc = pickle.load(open('dod5_xi36_closed.pkl','rb'))
Rc = pickle.load(open('dod5_R_content.pkl','rb'))
aco = pickle.load(open('dod5_atom_coefs_all.pkl','rb'))
sd = pickle.load(open(D1+"S_atom_decomp.pkl","rb"))
atoms = pickle.load(open(D1+"atoms620.pkl","rb"))
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
IDX=[1,2,5,6]
from fractions import Fraction as Fr
def fadd(a,b): return (a[0]+b[0], a[1]+b[1])
def fmul(a,b): return (a[0]*b[0]+3*a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def fscale(a,q): return (a[0]*q, a[1]*q)
ZERO=(Fr(0),Fr(0))
def pair_of_pq(pq):
    p,q = pq
    return (Fr(p.numerator, p.denominator), Fr(q.numerator, q.denominator))
gam={}; dl={}
for i in IDX:
    for j in IDX:
        if i==j: continue
        gam[(i,j)] = pdiv(Xs[j], psub(Xs[j], Xs[i]))
        dl[(i,j)]  = pmul(gam[(i,j)], psub(ONE, Xs[i]))

# ---- raw Li3 coefficients in S: dict pt -> pair ----
rawc = {}
def acc_pt(pt, c):
    rawc[pt] = fadd(rawc.get(pt, ZERO), c)
# Xi36: A3 at basis args
for j,z in enumerate(BASIS3_ZS):
    if xc['A3'][j] != ZERO: acc_pt(z, xc['A3'][j])
# R: content coeffs
for pt,c in Rc['li3_content'].items():
    acc_pt(pt, pair_of_pq(psplit(c)))
# L3x, Ax
for i in IDX:
    c = pair_of_pq(aco[('L3x',i)]['1'])
    acc_pt(Xs[i], c)
    ca = pair_of_pq(aco[('Ax',i)]['1'])
    acc_pt(Xs[i], ca)
    acc_pt(psub(ONE,Xs[i]), (Fr(-ca[0]), Fr(-ca[1])))
# T3p: c*(2/x_i) at 1-dl, c*(-2/x_i) at 1-gam
for i in IDX:
    for j in IDX:
        if i==j: continue
        c = pair_of_pq(aco[('T3p',(i,j))]['1'])
        two_over_xi = psplit(pdiv(P(2,0), Xs[i]))
        c2 = fmul(c, two_over_xi)
        acc_pt(psub(ONE, dl[(i,j)]), c2)
        acc_pt(psub(ONE, gam[(i,j)]), (Fr(-c2[0]), Fr(-c2[1])))
print("raw Li3 points in S:", len(rawc))
for pt in sorted(rawc, key=pfloat):
    print("   %-14s ~%+.8f  (%s) + (%s) r3"%(str(pt), pfloat(pt), rawc[pt][0], rawc[pt][1]))
pickle.dump(rawc, open('dod5_rawli3_coeffs.pkl','wb'))
print("t=%.0fs"%(time.time()-t0))
