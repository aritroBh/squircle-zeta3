import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
t0=time.time()
mp.mp.dps = 620
cl = pickle.load(open('dod5_closures36.pkl','rb'))
xa = pickle.load(open('dod5_xi36_acc.pkl','rb'))
ac = pickle.load(open(D2+"atom_coefs_pq.pkl","rb"))
atoms = pickle.load(open(D1+"atoms620.pkl","rb"))
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}

# F-elements as (Fr, Fr) pairs: c = c0 + c1*sqrt3
def fadd(a,b): return (a[0]+b[0], a[1]+b[1])
def fmul(a,b): return (a[0]*b[0]+3*a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def fscale(a,q): return (a[0]*q, a[1]*q)
ZERO=(Fr(0),Fr(0))

# coefficient of each of the 36 terms
coefF = {}
for (kind, ij),(p,q) in ac.items():
    c = (Fr(p), Fr(q))
    if kind=='Lam':
        coefF[('Lam',ij)] = c
    else:
        i,j,k = ij
        xi_p, xi_q = psplit(Xs[i])
        den = xi_p*xi_p - 3*xi_q*xi_q
        gp = (Fr(p)*xi_p - 3*Fr(q)*xi_q)/den; gq = (Fr(q)*xi_p - Fr(p)*xi_q)/den
        coefF[('SwG',ij)] = (gp,gq)
        coefF[('SwD',ij)] = (-gp,-gq)

# accumulate: Li3[8] pairs, prod[19] pairs, invisible[5] pairs (F-basis {z3,L3c,z2l2,z2l3,z2lu})
A3 = [ZERO]*8; AP = [ZERO]*19; AI = [ZERO]*5
for key, clo in cl.items():
    c = coefF[key]
    for j in range(8):
        if clo['solT'][j]: A3[j] = fadd(A3[j], fscale(c, Fr(clo['solT'][j])))
    for j in range(19):
        if clo['solP'][j]: AP[j] = fadd(AP[j], fscale(c, Fr(clo['solP'][j])))
    rel = clo['rel']
    for i in range(10):
        qi = -Fr(rel[i+1])/Fr(rel[0])
        m, par = i//2, i%2   # basis10: even -> X_m, odd -> sqrt3 X_m
        add = (qi,Fr(0)) if par==0 else (Fr(0),qi)
        AI[m] = fadd(AI[m], fmul(c, add))

# numeric value of the closed form
C = consts(620)
LI3B = [Li3v(z) for z in BASIS3_ZS]
PBV = [prodval(m, C) for m in PROD_META]
IBV = [C['z3'], C['L3c'], C['z2']*C['l2'], C['z2']*C['l3'], C['z2']*C['lu']]
def fval(p):
    return mp.mpf(p[0].numerator)/p[0].denominator + mp.sqrt(3)*mp.mpf(p[1].numerator)/p[1].denominator
tot = mp.mpf(0)
for j in range(8): tot += fval(A3[j])*LI3B[j]
for j in range(19): tot += fval(AP[j])*PBV[j]
for m in range(5): tot += fval(AI[m])*IBV[m]
target = xa['Xi36_acc']
print("closed-form Xi36 =", mp.nstr(tot, 40))
print("accurate Xi36    =", mp.nstr(target, 40))
print("|diff| =", mp.nstr(abs(tot-target),3))
print()
print("Li3 coefficients (8 basis args):")
for j in range(8):
    if A3[j]!=ZERO: print("   Li3(%-12s ~%+.8f): %s + (%s)*sqrt3"%(str(BASIS3_ZS[j]), pfloat(BASIS3_ZS[j]), A3[j][0], A3[j][1]))
print("product coefficients (19):")
for j in range(19):
    if AP[j]!=ZERO: print("   prod[%s]: %s + (%s)*sqrt3"%(PROD_META[j], AP[j][0], AP[j][1]))
print("invisible coefficients (F-basis z3, L3chi, z2l2, z2l3, z2lu):")
nms = ['z3','L3chi','z2*l2','z2*l3','z2*lu']
for m in range(5):
    if AI[m]!=ZERO: print("   %s: %s + (%s)*sqrt3"%(nms[m], AI[m][0], AI[m][1]))
pickle.dump({'A3':A3,'AP':AP,'AI':AI,'val':tot}, open('dod5_xi36_closed.pkl','wb'))
print("saved dod5_xi36_closed.pkl, elapsed %.0fs"%(time.time()-t0))
