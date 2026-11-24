
script = r'''
# step_11.py -- numeric function identity: D(c) = Omega - DT_A - Classical at c=1/2, 1/3
import sys, math, pickle, re
sys.path.insert(0,'.')
import mpmath as mp
import sympy as sp
from fractions import Fraction as Fr
mp.mp.dps = 50
c_sym = sp.symbols('c')

def hpl_coeffs(w, N):
    if len(w) == 0:
        cc = [mp.mpf(0)]*(N+1); cc[0]=mp.mpf(1); return cc
    a, rest = w[0], w[1:]
    cc = hpl_coeffs(rest, N)
    out = [mp.mpf(0)]*(N+1)
    if a == 0:
        for n in range(1,N+1): out[n] = cc[n]/n
    elif a == 1:
        S = mp.mpf(0)
        for n in range(1,N+1):
            S += cc[n-1]; out[n] = S/n
    else:
        A = mp.mpf(0)
        for n in range(1,N+1):
            A = cc[n-1]-A; out[n] = A/n
    return out
def H(w, x):
    w = tuple(w)
    if w == (0,): return mp.log(x)
    N = int((mp.mp.dps+15)/abs(math.log10(abs(x)))) + 10
    cc = hpl_coeffs(w, N)
    tot = mp.mpf(0); xn = x
    for n in range(1,N+1):
        tot += cc[n]*xn; xn *= x
    return tot
def Li22s(x,y):
    tot = mp.mpf(0); inner = mp.mpf(0); xn = x; n = 1
    while True:
        if n>=2:
            m = n-1; inner += y**m/mp.mpf(m)**2
        term = xn*inner/mp.mpf(n)**2
        tot += term
        if n>40 and abs(term)*(1+x/(1-x)*2) < mp.mpf(10)**(-mp.mp.dps-5): break
        n += 1; xn *= x
    return tot

# Omega from HPL list
import json
dd = json.load(open('../../10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
def parse(term):
    m = re.match(r'c\^(-?\d+)\s*(.*)', term[0]); e = int(m.group(1)); rest = m.group(2).strip()
    w = [int(x) for x in rest[2:-1].split(',')] if rest.startswith('H') else None
    return e, w
def Omega(c0):
    tot = mp.mpf(0)
    for t in dd:
        e,w = parse(t)
        if w is None or len(w)!=4: continue
        tot += (mp.mpf(c0)**e) * mp.mpf(Fr(t[1])) * H(w, c0)
    return tot
def DTA(c0):
    return ((1-c0)/c0)*Li22s(mp.mpf('0.5'),1-c0) - ((1+c0)/c0)*Li22s(mp.mpf('0.5'),1+c0)

# classical combination
sols2 = pickle.load(open('classical_sols_tracked.pkl','rb'))
seen = pickle.load(open('li_arg_pairs.pkl','rb'))
gx = {}; g1 = {}
for (u,x), arg in seen.items():
    gx.setdefault((u,x), sp.sympify(arg))
for (u,x), arg in seen.items():
    g1.setdefault((u,), sp.sympify(arg))
def ev_Li(n, arg, c0):
    z = sp.N(arg.subs(c_sym, sp.Rational(c0)), 30)
    return mp.polylog(n, mp.mpc(z))
def ev_Li1(arg, c0):
    z = sp.N(arg.subs(c_sym, sp.Rational(c0)), 30)
    return -mp.log(1-mp.mpc(z))
def ev_name(name, c0):
    kind = name[0]
    if kind=='Li4': return ev_Li(4, gx[name[1]], c0)
    if kind=='L3L1': return ev_Li(3, gx[name[1]], c0)*ev_Li1(g1[name[2]], c0)
    if kind=='L2L2': return ev_Li(2, gx[name[1]], c0)*ev_Li(2, gx[name[2]], c0)
    if kind=='L2L1L1': return ev_Li(2, gx[name[1]], c0)*ev_Li1(g1[name[2]], c0)*ev_Li1(g1[name[3]], c0)
    if kind=='L1^4':
        p = mp.mpf(1)
        for k in name[1]: p *= ev_Li1(g1[k], c0)
        return p
def Classical(c0):
    tot = mp.mpc(0)
    for e in (-1,0,1):
        for name, cf in sols2[e].items():
            tot += (mp.mpf(c0)**e) * mp.mpf(cf) * ev_name(name, c0)
    return tot

for c0s in ['1/2','1/3']:
    c0 = mp.mpf(c0s)
    Om = Omega(c0); D = DTA(c0); Cl = Classical(c0)
    resid = Om - D - Cl
    print(f"--- c = {c0s}")
    print("  Omega   =", mp.nstr(Om, 45))
    print("  DT_A    =", mp.nstr(D, 45))
    print("  Class   =", mp.nstr(Cl, 45))
    print("  residual= ", mp.nstr(resid, 30))
'''
open('steps/step_11.py','w').write(script)
import subprocess, time
t=time.time()
r = subprocess.run(['nice','-n','15','<project>/.venv/bin/python','steps/step_11.py'],capture_output=True,text=True,timeout=1400)
print(r.stdout[-4000:]); print(r.stderr[-3000:]); print("elapsed", time.time()-t)
open('progress.log','a').write("step_11 first numeric identity residual at c=1/2,1/3\n")
