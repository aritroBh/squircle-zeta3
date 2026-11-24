import sympy as sp, pickle, os, math, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
r5 = sp.sqrt(5)
g = pickle.load(open('geom_exact.pkl','rb'))
apol, appol, h2oL2 = g['apol'], g['appol'], g['h2oL2']
Z = sp.Integer(0)
def prune(P):
    return {k:v for k,v in P.items() if any(x != 0 for x in v)}
def eadd(P,Q):
    R = dict(P)
    for k in Q:
        if k in R: R[k] = [sp.expand(R[k][i]+Q[k][i]) for i in range(4)]
        else: R[k] = list(Q[k])
    return R
def ekmul(P,Q):
    R = {}
    for k1,v1 in P.items():
        s1 = sum(k1)
        for k2,v2 in Q.items():
            k = tuple(k1[i]+k2[i] for i in range(5))
            if s1+sum(k2) > 3: continue
            c = [Z]*4
            for i in range(4):
                for j in range(4-i):
                    if v1[i] != 0 and v2[j] != 0:
                        c[i+j] = sp.expand(c[i+j]+v1[i]*v2[j])
            if k in R: R[k] = [sp.expand(R[k][i]+c[i]) for i in range(4)]
            else: R[k] = c
    return R
def escale(P,s):
    return {k:[sp.expand(v[i]*s) for i in range(4)] for k,v in P.items()}
ONE = {(0,0,0,0,0):[sp.Integer(1),Z,Z,Z]}
def pmul(P,Q):
    R = [{} for _ in range(4)]
    for i in range(4):
        for j in range(4-i):
            R[i+j] = eadd(R[i+j], ekmul(P[i],Q[j]))
    return R
def pdiv(P,Q):
    R = [{} for _ in range(4)]
    q0 = Q[0][(0,0,0,0,0)][0]
    for i in range(4):
        acc = P[i]
        for j in range(i):
            acc = eadd(acc, escale(ekmul(R[j],Q[i-j]), sp.Integer(-1)))
        R[i] = {k:[sp.expand(x/q0) for x in v] for k,v in acc.items()}
    return R
def plog(P):
    Pp = [escale(P[i], sp.Integer(i)) for i in range(1,4)] + [{}]
    Qd = pdiv(Pp,P)
    L = [{}]
    for i in range(3): L.append(escale(Qd[i], sp.Rational(1,i+1)))
    return L
def pexp(P):
    E = [ONE] + [{} for _ in range(3)]
    for k in range(1,4):
        acc = {}
        for j in range(1,k+1):
            acc = eadd(acc, escale(ekmul(P[j],E[k-j]), sp.Rational(j)))
        E[k] = escale(acc, sp.Rational(1,k))
    return E
t0=time.time()
F = [{} for _ in range(4)]; S = [{} for _ in range(4)]
for k in range(4):
    fc = sp.Rational(1, math.factorial(k))
    for i in range(5):
        key = tuple(k if j==i else 0 for j in range(5))
        F[k] = eadd(F[k], {key:[sp.expand(apol[i][0]*fc), sp.expand(apol[i][1]*fc), Z, Z]})
        S[k] = eadd(S[k], {key:[sp.expand(appol[i][0]*fc), sp.expand(appol[i][1]*fc), Z, Z]})
F = [prune(x) for x in F]; S = [prune(x) for x in S]
assert F[0] == ONE, F[0]
assert S[0] == {}, S[0]
print("F0=1, S0=0  [EXACT]", time.time()-t0)
Phi = pdiv(S,F); lg = plog(F)
inv = [escale(ONE, sp.Integer(s)) for s in [1,-1,1,-1]]
m2lg = [escale(lg[i], sp.Integer(-2)) for i in range(4)]
facE = pexp(pmul(m2lg, inv))
Phisq = pmul(Phi,Phi)
o2mt = {(0,0,0,0,0):[sp.Integer(1), sp.Integer(-2), Z, Z]}
yypE = {(0,0,0,0,0):[sp.expand(h2oL2+sp.Rational(1,4)), sp.Integer(-1), sp.Integer(1), Z]}
Q = [{} for _ in range(4)]
Q[0] = ONE
for k in range(1,4): Q[k] = ekmul(Phi[k], o2mt)
for k in range(4):
    Q[k] = eadd(Q[k], ekmul(Phisq[k], yypE))
G = pmul(facE, Q)
sq = [{} for _ in range(4)]; sq[0] = ONE
for k in range(1,4):
    acc = G[k]
    for j in range(1,k):
        acc = eadd(acc, escale(ekmul(sq[j],sq[k-j]), sp.Integer(-1)))
    sq[k] = escale(acc, sp.Rational(1,2))
E3n = prune(sq[3])
print("engine done", time.time()-t0, "monomials:", len(E3n))
def pq(e):
    e = sp.expand(e)
    p = sp.simplify((e + e.subs(r5,-r5))/2)
    q = sp.simplify((e - e.subs(r5,-r5))/(2*r5))
    return (sp.Rational(p), sp.Rational(q))
E3n_pq = {k:[pq(x) for x in v] for k,v in E3n.items()}
pickle.dump(E3n_pq, open('E3n_sym_pq.pkl','wb'))
print("saved; sample:")
for k in list(E3n_pq)[:4]:
    print(k, E3n_pq[k])
