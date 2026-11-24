import mpmath as mp, sympy as sp, pickle
mp.mp.dps = 50
N = 256; R = mp.mpf('0.35')
om = [mp.exp(2*mp.pi*mp.j*k/N) for k in range(N)]
def Ek_cauchy(L, h, s0, t, kmax=5):
    l0 = mp.log(t); l1 = mp.log(1-t); sig = s0 + L*t
    vals = []; aF = aw = None
    for j in range(N):
        d = R*om[j]; p = 1 + d
        e0 = mp.exp(d*l0); e1 = mp.exp(d*l1)
        F = t*e0 + (1-t)*e1
        Phi = (e0 - e1)/F
        w2 = (L - sig*Phi)**2 + h**2*Phi**2
        def clog(z, pa):
            a = mp.arg(z)
            if pa is not None: a += 2*mp.pi*int(mp.nint((pa - a)/(2*mp.pi)))
            return mp.log(abs(z)) + 1j*a, a
        logF, aF = clog(F, aF); logw, aw = clog(w2, aw)
        vals.append(mp.exp(logw/2 - logF/p))
    return [mp.re(mp.fsum(vals[j]*om[j]**(-k) for j in range(N))/(N*R**k)) for k in range(kmax+1)]

# symbolic reference values of E_k(t) at several t, edge A
tt, al, be, Ls, hs, s0s, l0s, l1s = sp.symbols('t alpha beta L h s0 l0 l1')
Rp = pickle.load(open('steps/step_01_reduced.pkl','rb'))
loc = {'t': tt, 'alpha': al, 'beta': be, 'L': Ls, 'h': hs, 's0': s0s, 'l0': l0s, 'l1': l1s}
E = [sp.sympify(s, locals=loc) for s in Rp['E']]
L = mp.sqrt(5); h = 4/mp.sqrt(5); s0 = -2/mp.sqrt(5)
import time
worst = mp.mpf(0)
for tq in [(1,100),(1,10),(3,10),(1,2),(7,10),(9,10),(99,100)]:
    t = mp.mpf(tq[0])/tq[1]
    Ec = Ek_cauchy(L, h, s0, t)
    sub = {al:1, be:1, Ls: sp.sqrt(5), hs: 4/sp.sqrt(5), s0s: -2/sp.sqrt(5),
           tt: sp.Rational(*tq), l0s: sp.log(sp.Rational(*tq)), l1s: sp.log(1-sp.Rational(*tq))}
    Esym = [mp.mpf(sp.N(c.subs(sub), 50)) for c in E]
    dmax = max(abs(Ec[k]-Esym[k]) for k in range(6))
    worst = max(worst, dmax)
    print(f"t={mp.nstr(t,4)}: max |Cauchy - sym| over k=0..5: {mp.nstr(dmax,3)}")
print("overall worst:", mp.nstr(worst,3))
