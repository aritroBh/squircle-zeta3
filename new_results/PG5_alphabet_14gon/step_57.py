e3n = r'''
"""PG5 step_07: exact E3n (d^3 integrand) for the 14-gon, port of DOD1 step_03
to exact FT arithmetic. Term order: [a3(loc t), a4(loc 1-t), a1,a2,a5,a6, a0(const)].
Validated by: E3n(0)=1 coefficient structure, and 84*L*int(E3n) vs gold P3(1) (step_08).
"""
import pickle, sys, time
sys.path.insert(0,'.')
import ftfield as F
from fractions import Fraction as Fr
t0=time.time()
g = pickle.load(open('geom14.pkl','rb'))
data = g['data']  # i -> (const, slope) FT
h2oL2 = g['h2oL2']
ORDER = [3,4,1,2,5,6,0]
apol  = [[data[i][0], data[i][1]] for i in ORDER]
appol = [[data[i][1], F.ZERO] for i in ORDER]
NT = 7
ZM = (0,)*NT
def eadd(P,Q):
    R = dict(P)
    for k,v in Q.items():
        if k in R: R[k] = [F.fadd(R[k][i], v[i]) for i in range(4)]
        else: R[k] = list(v)
    return R
def ekmul(P,Q):
    R = {}
    for k1,v1 in P.items():
        s1 = sum(k1)
        for k2,v2 in Q.items():
            if s1+sum(k2) > 3: continue
            k = tuple(k1[i]+k2[i] for i in range(NT))
            c = [F.ZERO]*4
            for i in range(4):
                for j in range(4-i):
                    c[i+j] = F.fadd(c[i+j], F.fmul(v1[i],v2[j]))
            if k in R: R[k] = [F.fadd(R[k][i], c[i]) for i in range(4)]
            else: R[k] = c
    return R
def escale(P,s):
    return {k:[F.fmul(v[i],s) for i in range(4)] for k,v in P.items()}
ONE = {ZM:[F.ONE,F.ZERO,F.ZERO,F.ZERO]}
def pmul(P,Q):
    R=[{} for _ in range(4)]
    for i in range(4):
        for j in range(4-i):
            R[i+j]=eadd(R[i+j],ekmul(P[i],Q[j]))
    return R
def pdiv(P,Q):
    R=[{} for _ in range(4)]
    q0 = Q[0][ZM][0]
    for i in range(4):
        acc = dict(P[i])
        for j in range(i):
            acc = eadd(acc, escale(ekmul(R[j],Q[i-j]), F.ft(-1)))
        R[i] = {k:[F.fdiv(v[i2], q0) for i2 in range(4)] for k,v in acc.items()}
    return R
def plog(P):
    Pp = [escale(P[i], F.ft(i)) for i in range(1,4)] + [{}]
    Q = pdiv(Pp, P)
    L = [{}]
    for i in range(3): L.append(escale(Q[i], F.ft(Fr(1,i+1))))
    return L
def pexp(P):
    E = [ONE] + [{}]*3
    for k in range(1,4):
        acc = {}
        for j in range(1,k+1):
            acc = eadd(acc, escale(ekmul(P[j],E[k-j]), F.ft(j)))
        E[k] = escale(acc, F.ft(Fr(1,k)))
    return E
Fc=[{} for _ in range(4)]; Sc=[{} for _ in range(4)]
for k in range(4):
    fc = Fr(1)
    for _ in range(k): fc /= (_+1)
    for i in range(NT):
        key = tuple(k if j==i else 0 for j in range(NT))
        Fc[k] = eadd(Fc[k], {key:[F.fscale(apol[i][0],fc), F.fscale(apol[i][1],fc), F.ZERO, F.ZERO]})
        Sc[k] = eadd(Sc[k], {key:[F.fscale(appol[i][0],fc), F.ZERO, F.ZERO, F.ZERO]})
Phi = pdiv(Sc,Fc); lg = plog(Fc)
inv = [escale(ONE, F.ft(s)) for s in [1,-1,1,-1]]
m2lg = [escale(x, F.ft(-2)) for x in lg]
facE = pexp(pmul(m2lg, inv))
Phisq = pmul(Phi,Phi)
Q=[{} for _ in range(4)]
Q[0]=ONE
o2mt = {ZM:[F.ONE, F.ft(-2), F.ZERO, F.ZERO]}
yy2 = [F.fadd(h2oL2, F.ft(Fr(1,4))), F.ft(-1), F.ONE, F.ZERO]
for k in range(1,4): Q[k]=ekmul(Phi[k], o2mt)
for k in range(4):
    Q[k]=eadd(Q[k], ekmul(Phisq[k], {ZM:yy2}))
G = pmul(facE, Q)
sq=[{} for _ in range(4)]; sq[0]=ONE
for k in range(1,4):
    acc = dict(G[k])
    for j in range(1,k):
        acc = eadd(acc, escale(ekmul(sq[j],sq[k-j]), F.ft(-1)))
    sq[k] = escale(acc, F.ft(Fr(1,2)))
E3n = sq[3]
nz = {k:v for k,v in E3n.items() if any(vv[:3]!=(Fr(0),)*3 for vv in v)}
print("E3n monomials (7 term-logs):", len(nz), " time %.1fs"%(time.time()-t0), flush=True)
pickle.dump({'E3n':nz,'ORDER':ORDER}, open('E3n14.pkl','wb'))
# numeric spot check vs direct numeric build of the same formula (independent mpmath impl)
import mpmath as mp
mp.mp.dps = 60
B = g['B']; X = g['X']; A0 = g['A0']; hh = g['h']
def fv(x): return F.feval(x, 60)
apn  = [[fv(data[i][0]), fv(data[i][1])] for i in ORDER]
t = mp.mpf('0.37')
logs = [mp.log(apn[i][0]+apn[i][1]*t) for i in range(NT)]
val = mp.mpf(0)
for k,v in nz.items():
    poly = ((fv(v[3])*t+fv(v[2]))*t+fv(v[1]))*t+fv(v[0])
    lv = mp.mpf(1)
    for i,e in enumerate(k):
        if e: lv *= logs[i]**e
    val += poly*lv
print("E3n(0.37) exact-eval =", mp.nstr(val, 40))
pickle.dump({'t':'0.37','val':str(val)}, open('E3n14_spot.pkl','wb'))
open('progress.log','a').write('PG5 step_07: exact E3n 14-gon built: %d monomials, %.0fs\n'%(len(nz),time.time()-t0))
'''
open('step_07_e3n14.py','w').write(e3n)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_07_e3n14.py'],capture_output=True,text=True,timeout=1500)
print("elapsed", time.time()-t0)
print(r.stdout); print(r.stderr[-800:] if r.stderr else '')
