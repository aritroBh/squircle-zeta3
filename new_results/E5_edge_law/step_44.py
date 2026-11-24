import mpmath as mp, sympy as sp, pickle
mp.mp.dps = 80
N = 192
om = [mp.exp(2*mp.pi*mp.j*k/N) for k in range(N)]
def cauchy(L, h, s0, t, R, kmax=5):
    l0 = mp.log(t); l1 = mp.log(1-t); sig = s0 + L*t
    vals = []; aF = aw = None; M = mp.mpf(0)
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
        g = mp.exp(logw/2 - logF/p); M = max(M, abs(g))
        vals.append(g)
    E = [mp.re(mp.fsum(vals[j]*om[j]**(-k) for j in range(N))/(N*R**k)) for k in range(kmax+1)]
    return E, M
def Ek2(t, L, h, s0):
    R = mp.mpf('0.35'); E1, M1 = cauchy(L, h, s0, t, R)
    for _ in range(6):
        E2, M2 = cauchy(L, h, s0, t, R/2)
        if max(abs(E1[k]-E2[k]) for k in range(6)) < mp.mpf('1e-64'): return E2, M2, R/2
        R = R/2; E1, M1 = E2, M2
    return E1, M1, R
tt, al, be, Ls, hs, s0s, l0s, l1s = sp.symbols('t alpha beta L h s0 l0 l1')
Rp = pickle.load(open('steps/step_01_reduced.pkl','rb'))
loc = {'t': tt, 'alpha': al, 'beta': be, 'L': Ls, 'h': hs, 's0': s0s, 'l0': l0s, 'l1': l1s}
E = [sp.sympify(s, locals=loc) for s in Rp['E']]
L = mp.sqrt(5); h = 4/mp.sqrt(5); s0 = -2/mp.sqrt(5)
for dt in ['1e-10', '1e-20', '1e-25']:
    for side in [0, 1]:
        t = mp.mpf(dt) if side == 0 else 1 - mp.mpf(dt)
        Ev, M, Rused = Ek2(t, L, h, s0)
        ts = sp.mpf(dt) if side == 0 else 1 - sp.mpf(dt)
        sub = {al:1, be:1, Ls: sp.sqrt(5), hs: 4/sp.sqrt(5), s0s: -2/sp.sqrt(5),
               tt: ts, l0s: sp.log(ts), l1s: sp.log(1-ts)}
        Esym = [mp.mpf(sp.N(c.subs(sub), 80)) for c in E]
        dmax = max(abs(Ev[k]-Esym[k]) for k in range(1, 6))
        print(f"t={'1-' if side else ''}{dt}: |E5|={mp.nstr(abs(Esym[5]),3)} M={mp.nstr(M,3)} R={mp.nstr(Rused,3)}"
              f" max err k=1..5: {mp.nstr(dmax,3)}")
