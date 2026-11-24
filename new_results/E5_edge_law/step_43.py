code = r'''# step_07b.py: convergence study of the per-t-Cauchy + tanh-sinh scheme on edge A only.
import mpmath as mp, time, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 80
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
    for _ in range(5):
        E2 = cauchy(L, h, s0, t, R/2)
        if max(abs(E1[k]-E2[k]) for k in range(6)) < mp.mpf('1e-64'): return E2
        R = R/2; E1 = E2
    return E1
def run_edge(L, h, s0, du0=mp.mpf('0.09'), umax=mp.mpf('3.6'), levels=6):
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
        print(f"  level {lev}: nodes<= {2*jmax+1}  E0-int err = {mp.nstr(abs(sums[0]-L),3)}  E5-int = {mp.nstr(sums[5],12)}")
    return sums
t0 = time.time()
V0, V1 = (mp.mpf(2), mp.mpf(0)), (mp.mpf(1), mp.mpf(2))
L = mp.sqrt((V1[0]-V0[0])**2+(V1[1]-V0[1])**2)
h = abs(V0[0]*V1[1]-V0[1]*V1[0])/L
tau = ((V1[0]-V0[0])/L, (V1[1]-V0[1])/L)
s0 = V0[0]*tau[0]+V0[1]*tau[1]
sums = run_edge(L, h, s0)
print("wall %.0f s" % (time.time()-t0))
'''
open('steps/step_07b.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_07b.py','80'],capture_output=True,text=True,timeout=5400)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-1500:])
