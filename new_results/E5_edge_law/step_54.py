code = r'''# step_10.py: INDEPENDENT numeric P^(5)(1) for the rhombus |x|^p+|y/2|^p<=1 (4 congruent edges),
# same machine as step_07 (per-t Cauchy + tanh-sinh, no symbolic E-series used), then
# PSLQ of P^(5)(1)/sqrt5 on {1, z2, z3, z2^2, z2z3, z5} (7 elements).
import mpmath as mp, time, json
mp.mp.dps = 100
N = 192
om = [mp.exp(2*mp.pi*mp.j*k/N) for k in range(N)]
def cauchy(L, h, s0, t, R, kmax=5):
    l0 = mp.log(t); l1 = mp.log(1-t); sig = s0 + L*t
    vals = []; aF = aw = None
    for j in range(N):
        d = R*om[j]; p = 1 + d
        e0 = mp.exp(d*l0); e1 = mp.exp(d*l1)
        F = t*e0 + (1-t)*e1
        Phi = (e0 - e1)/F
        w2 = (L - sig*Phi)**2 + h**2*Phi**2
        a1 = mp.arg(F); a2 = mp.arg(w2)
        if aF is not None:
            a1 += 2*mp.pi*int(mp.nint((aF - a1)/(2*mp.pi)))
            a2 += 2*mp.pi*int(mp.nint((aw - a2)/(2*mp.pi)))
        aF, aw = a1, a2
        logF = mp.log(abs(F)) + 1j*a1; logw = mp.log(abs(w2)) + 1j*a2
        vals.append(mp.exp(logw/2 - logF/p))
    return [mp.re(mp.fsum(vals[j]*om[j]**(-k) for j in range(N))/(N*R**k)) for k in range(kmax+1)]
def Ek(t, L, h, s0):
    R = mp.mpf('0.35'); E1 = cauchy(L, h, s0, t, R)
    for _ in range(6):
        E2 = cauchy(L, h, s0, t, R/2)
        if max(abs(E1[k]-E2[k]) for k in range(6)) < mp.mpf('1e-66'): return E2
        R = R/2; E1 = E2
    return E1
def tanh_sinh_edge(L, h, s0, du0=mp.mpf('0.09'), umax=mp.mpf('4.6'), levels=5):
    sums = None; prev = None
    for lev in range(levels):
        du = du0/2**lev
        jmax = int(mp.ceil(umax/du))
        js = range(-jmax, jmax+1) if lev == 0 else [j for j in range(-jmax, jmax+1) if j % 2]
        if sums is None: sums = [mp.mpf(0)]*6
        else: sums = [s/2 for s in sums]
        for j in js:
            u = j*du
            sh = mp.sinh(u); th = mp.tanh(mp.pi/2*sh)
            t = (1+th)/2
            w = du*(mp.pi/4)*mp.cosh(u)*(1/mp.cosh(mp.pi/2*sh))**2
            Ev = Ek(t, L, h, s0)
            for k in range(6): sums[k] += w*Ev[k]
        print(f"    level {lev}: E0 err={mp.nstr(abs(sums[0]-L),3)}")
    return sums
t0 = time.time()
V0, V1 = (mp.mpf(1), mp.mpf(0)), (mp.mpf(0), mp.mpf(2))
L = mp.sqrt((V1[0]-V0[0])**2+(V1[1]-V0[1])**2)
h = abs(V0[0]*V1[1]-V0[1]*V1[0])/L
tau = ((V1[0]-V0[0])/L, (V1[1]-V0[1])/L)
s0 = V0[0]*tau[0]+V0[1]*tau[1]
sums = tanh_sinh_edge(L, h, s0)
Pk = [mp.factorial(k)*4*sums[k] for k in range(6)]
z2, z3, z5 = mp.zeta(2), mp.zeta(3), mp.zeta(5)
r5 = mp.sqrt(5)
exact5 = r5*(1152*z5/25 - 768*z2*z3/25 - 4608*z2**2/125 + 2304*z3/5 - 1536*z2/5)
print("wall %.0f s" % (time.time()-t0))
print("P^(5)(1) numeric =", mp.nstr(Pk[5], 60))
print("P^(5)(1) exact   =", mp.nstr(exact5, 60))
print("|diff| =", mp.nstr(abs(Pk[5]-exact5), 3))
json.dump({'P5': mp.nstr(Pk[5], 75)}, open('steps/step_10_num.json', 'w'))
# PSLQ: V = P5/sqrt5 on {1, z2, z3, z2^2, z2z3, z5}
V = Pk[5]/r5
basis = [mp.mpf(1), z2, z3, z2**2, z2*z3, z5]
bn = ['1', 'z2', 'z3', 'z2^2', 'z2z3', 'z5']
rel = mp.pslq(mp.matrix([V] + basis), tol=mp.mpf('1e-60'), maxcoeff=10**7, maxsteps=10000)
print("PSLQ relation:", rel)
if rel is not None and rel[0] != 0:
    q = [mp.mpf(rel[i+1])/(-rel[0]) for i in range(6)]
    print("V =", {bn[i]: str(sp_q) for i, sp_q in enumerate([mp.nstr(qq, 10) for qq in q]) if qq != 0})
    print("check exact coeffs [0, -1536/5, 2304/5, -4608/125, -768/25, 1152/25]")
    from fractions import Fraction
    print("as rationals:", [str(Fraction(int(mp.mpf(rel[i+1])), int(-rel[0]))) for i in range(6)])
'''
open('steps/step_10.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_10.py'],capture_output=True,text=True,timeout=7200)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-2000:])
