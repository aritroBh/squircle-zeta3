import mpmath as mp, sympy as sp, pickle
mp.mp.dps = 50
N = 220; R = mp.mpf('0.5')
om = [mp.exp(2*mp.pi*mp.j*k/N) for k in range(N)]

def Ek_cauchy(L, h, s0, t, kmax=5, track=True):
    l0 = mp.log(t); l1 = mp.log(1-t); sig = s0 + L*t
    vals = []; prev = None
    for j in range(N):
        d = R*om[j]; p = 1 + d
        e0 = mp.exp(d*l0); e1 = mp.exp(d*l1)
        F = t*e0 + (1-t)*e1
        Phi = (e0 - e1)/F
        w2 = (L - sig*Phi)**2 + h**2*Phi**2
        if not track:
            vals.append(mp.sqrt(w2)*mp.exp(-mp.log(F)/p)); continue
        # continuous branch: unwind arg of F and w2 relative to previous node
        def clog(z, prev_arg):
            a = mp.arg(z)
            if prev_arg is not None:
                a += 2*mp.pi*round((prev_arg - a)/(2*mp.pi))
            return mp.log(abs(z)) + 1j*a, a
        if prev is None:
            logF, aF = clog(F, None); logw, aw = clog(w2, None)
        else:
            logF, aF = clog(F, prev[0]); logw, aw = clog(w2, prev[1])
        prev = (aF, aw)
        vals.append(mp.exp(logw/2 - logF/p))
    return [mp.re(mp.fsum(vals[j]*om[j]**(-k) for j in range(N))/(N*R**k)) for k in range(kmax+1)]

# edge A geometry
L = mp.sqrt(5); h = 4/mp.sqrt(5); s0 = -2/mp.sqrt(5)
for tv in ['0.3', '0.7']:
    t = mp.mpf(tv)
    E_plain = Ek_cauchy(L, h, s0, t, track=False)
    E_track = Ek_cauchy(L, h, s0, t, track=True)
    # symbolic reference from step_01 E-series
    tt, al, be, Ls, hs, s0s, l0s, l1s = sp.symbols('t alpha beta L h s0 l0 l1')
    Rp = pickle.load(open('steps/step_01_reduced.pkl','rb'))
    loc = {'t': tt, 'alpha': al, 'beta': be, 'L': Ls, 'h': hs, 's0': s0s, 'l0': l0s, 'l1': l1s}
    E = [sp.sympify(s, locals=loc) for s in Rp['E']]
    sub = {al:1, be:1, L:sp.sqrt(5), h:4/sp.sqrt(5), s0:-2/sp.sqrt(5), tt:sp.Rational(3,10) if tv=='0.3' else sp.Rational(7,10),
           l0s: sp.log(sp.Rational(3,10) if tv=='0.3' else sp.Rational(7,10)),
           l1s: sp.log(sp.Rational(7,10) if tv=='0.3' else sp.Rational(3,10))}
    Esym = [sp.N(c.subs(sub), 50) for c in E]
    print("t =", tv)
    for k in range(6):
        print(f"  E{k}: plain={mp.nstr(E_plain[k],12)} track={mp.nstr(E_track[k],12)} sym={mp.nstr(Esym[k],12)}")
