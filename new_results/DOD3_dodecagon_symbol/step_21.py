import sys, os, pickle, time
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from dod3_symbol import *

base = os.path.dirname(os.getcwd())
ac = pickle.load(open(os.path.join(base,'DOD2_dodecagon_finish','atom_coefs_pq.pkl'),'rb'))
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
IDX=[1,2,5,6]
gam = {}; dl = {}
for i in IDX:
    for j in IDX:
        if i==j: continue
        gam[(i,j)] = pdiv(Xs[j], psub(Xs[j], Xs[i]))
        dl[(i,j)]  = pmul(gam[(i,j)], psub(ONE, Xs[i]))

# ---- assemble the 36 Li_{2,1} terms of Xi36 with exact F-coefficients ----
terms = []   # (kind, key, (x,y) args as pairs, coeff as pair in F)
for (kind, ij), (p,q) in ac.items():
    if kind=='Lam':
        i,j = ij
        c = P(p.numerator, q.numerator, p.denominator)  # p+q sqrt3 ; p,q have same denom here
        # careful: build exactly from fractions
        c = pscale(ONE, Fr(p))  # placeholder replaced below
        c = None
        m = P(p.numerator, q.numerator, 1)
        c = P(0,0)
        # c = p + q*r3 exactly:
        c = (Fr(p), Fr(q))
        xarg = Xs[i]; yarg = pdiv(Xs[j], Xs[i])
        terms.append(('Lam', (i,j), xarg, yarg, c))
    else:
        i,j,k = ij
        cpq = (Fr(p), Fr(q))
        xi = Xs[i]
        # coefficient for SwG = c/xi ; for SwD = -c/xi
        # c/xi in F: (p+q r3)/(xi) ; compute as pair then split
        xi_p, xi_q = psplit(xi)
        # (p+q r3)/(xi_p + xi_q r3) = (p+q r3)(xi_p - xi_q r3)/(xi_p^2-3 xi_q^2)
        den = xi_p*xi_p - 3*xi_q*xi_q
        gp = (p*xi_p - 3*q*xi_q)/den
        gq = (q*xi_p - p*xi_q)/den
        xg = pdiv(gam[(i,k)], gam[(i,j)]); yg = gam[(i,j)]
        xd = pdiv(dl[(i,k)], dl[(i,j)]);  yd = dl[(i,j)]
        terms.append(('SwG', (i,j,k), xg, yg, (gp, gq)))
        terms.append(('SwD', (i,j,k), xd, yd, (-gp, -gq)))

print("assembled %d Li21 terms"%len(terms))
# numeric cross-check of Xi36 vs stored: Xi = Xi36 + cT3a*R ; we can at least check the symbol later.
# ---- compute symbols S0 (1-part), S1 (sqrt3-part) ----
S0 = {}; S1 = {}
def acc(S, symb, a):
    for t,c in symb.items():
        v = S.get(t, Fr(0)) + a*c
        if v==0: S.pop(t,None)
        else: S[t]=v
for kind,key,xarg,yarg,(c0,c1) in terms:
    s = S_Li21(xarg, yarg)
    if c0: acc(S0, s, c0)
    if c1: acc(S1, s, c1)
print("S(Xi36_0): %d tensor terms;  S(Xi36_1): %d tensor terms"%(len(S0), len(S1)))
pickle.dump({'terms':[(k,kk,x,y,c) for k,kk,x,y,c in terms],'S0':S0,'S1':S1}, open('xi36_symbol.pkl','wb'))
print("saved xi36_symbol.pkl")

# ---- solve against Li3 gens + Decomp ----
dat = pickle.load(open('li3_gens.pkl','rb'))
zs, gens = dat['zs'], dat['gens']
Q = Quotient()
t0=time.time()
t0v = Q.coord(S0); t1v = Q.coord(S1)
ok0, c0s = Q.solve(t0v, gens)
ok1, c1s = Q.solve(t1v, gens)
print("solve Xi36_0:", ok0, " Xi36_1:", ok1, " (t=%.1fs)"%(time.time()-t0))
if ok0 and ok1:
    used = sorted(set([i for i,c in enumerate(c0s) if c!=0] + [i for i,c in enumerate(c1s) if c!=0]))
    print("number of Li3 arguments used:", len(used))
    print("\n%-4s %-16s %-22s %-22s"%("i","z","a0 (1-part)","a1 (sqrt3-part)"))
    for i in used:
        z = zs[i]
        print("%-4d %-16s %-22s %-22s"%(i, str(z), str(c0s[i]), str(c1s[i])))
    pickle.dump({'zs':zs,'c0':c0s,'c1':c1s,'used':used}, open('xi36_li3_solution.pkl','wb'))
    print("saved xi36_li3_solution.pkl")
