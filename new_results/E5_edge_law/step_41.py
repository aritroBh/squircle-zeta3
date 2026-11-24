code = r'''# step_07.py (v2): INDEPENDENT high-precision P^(k)(1), k=2..5, sheared parallelogram +-(2,0), +-(1,2).
# Per-t Cauchy extraction of E_k(t) = [d^k] |x_{1+d}'(t)| with continuous-branch logs and
# ADAPTIVE radius (halved until two radii agree): handles the w2-zeros approaching d=0 at the
# endpoints.  Manual tanh-sinh over t with level refinement; E_0-integral = L as built-in check.
import mpmath as mp, time, json, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 80
N = 256
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

def Ek(t, L, h, s0, tol):
    R = mp.mpf('0.35')
    E1 = cauchy(L, h, s0, t, R)
    for _ in range(6):
        E2 = cauchy(L, h, s0, t, R/2)
        if max(abs(E1[k]-E2[k]) for k in range(6)) < tol:
            return E2
        R = R/2; E1 = E2
    return E1

def tanh_sinh_edge(L, h, s0, tol, du0=mp.mpf('0.09'), umax=mp.mpf('3.6'), levels=5):
    prev = None; sums = None
    for lev in range(levels):
        du = du0/2**lev
        jmax = int(mp.ceil(umax/du))
        js = range(-jmax, jmax+1) if lev == 0 else [j for j in range(-jmax, jmax+1) if j % 2]
        if sums is None: sums = [mp.mpf(0)]*6
        for j in js:
            u = j*du
            sh = mp.sinh(u); th = mp.tanh(mp.pi/2*sh)
            t = (1+th)/2
            w = du*(mp.pi/4)*mp.cosh(u)*(1/mp.cosh(mp.pi/2*sh))**2
            Ev = Ek(t, L, h, s0, tol/100)
            for k in range(6): sums[k] += w*Ev[k]
        if prev is not None:
            diff = max(abs(sums[k]-prev[k]) for k in range(6))
            if diff < tol: break
        prev = list(sums)
    return sums, lev

t0 = time.time()
edges = [((2,0),(1,2)), ((1,2),(-2,0))]
tot = [mp.mpf(0)]*6
for (V0, V1) in edges:
    V0 = tuple(map(mp.mpf, V0)); V1 = tuple(map(mp.mpf, V1))
    L = mp.sqrt((V1[0]-V0[0])**2 + (V1[1]-V0[1])**2)
    h = abs(V0[0]*V1[1]-V0[1]*V1[0])/L
    tau = ((V1[0]-V0[0])/L, (V1[1]-V0[1])/L)
    s0 = V0[0]*tau[0]+V0[1]*tau[1]
    sums, lev = tanh_sinh_edge(L, h, s0, mp.mpf('1e-62'))
    print(f"edge V0={V0}: L={mp.nstr(L,6)} levels used={lev} E0-int={mp.nstr(sums[0],12)} (L={mp.nstr(L,12)})"
          f" E0 err={mp.nstr(abs(sums[0]-L),3)}  wall={time.time()-t0:.0f}s")
    for k in range(6): tot[k] += 2*sums[k]
z2, z3, z5 = mp.zeta(2), mp.zeta(3), mp.zeta(5)
r5, r13 = mp.sqrt(5), mp.sqrt(13)
exact = {
 2: z2*(64*r13/169 + 64*r5/25),
 3: z2*(-384*r5/25 - 384*r13/169) + z3*(192*r13/169 + 192*r5/25),
 4: z2**2*(-29952*r5/625 - 315648*r13/142805) + z2*(2304*r13/169 + 2304*r5/25) + z3*(-2304*r5/25 - 2304*r13/169),
 5: z2**2*(1262592*r13/28561 + 119808*r5/125) + z2*z3*(-10752*r5/25 - 821760*r13/28561)
    + z2*(-3072*r5/5 - 15360*r13/169) + z3*(23040*r13/169 + 4608*r5/5) + z5*(-6912*r5/25 - 172800*r13/28561)}
print("dps", mp.mp.dps, "wall %.0f s" % (time.time()-t0))
out = {}
for k in [2,3,4,5]:
    Pk = mp.factorial(k)*tot[k]
    out[str(k)] = mp.nstr(Pk, 70)
    print(f"P^({k})(1) = {mp.nstr(Pk, 40)}")
    print(f"   exact = {mp.nstr(exact[k], 40)}   |diff| = {mp.nstr(abs(Pk-exact[k]),3)}")
json.dump(out, open('steps/step_07_num.json', 'w'))
'''
open('steps/step_07.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_07.py','80'],capture_output=True,text=True,timeout=5400)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-2000:])
