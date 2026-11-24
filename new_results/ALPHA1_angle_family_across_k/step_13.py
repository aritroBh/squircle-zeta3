"""step_13.py -- PSLQ T(alpha), careful labels."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 70

def T_alpha(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

cases = {
  '1/2': (mp.mpf(1)/2, [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log3^2', mp.log(3)**2), ('log23', mp.log(2)*mp.log(3))],
          [(str(a), a) for a in [mp.mpf(1)/2, mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/4, mp.mpf(3)/4, -mp.mpf(1)/2, -mp.mpf(1)/3, -mp.mpf(2)/3, -mp.mpf(1)/4, mp.mpf(1)/9]]),
  '3/5': (mp.mpf(3)/5, [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log5^2', mp.log(5)**2), ('log25', mp.log(2)*mp.log(5))],
          [(str(a), a) for a in [mp.mpf(1)/2, mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/5, mp.mpf(2)/5, mp.mpf(3)/5, mp.mpf(4)/5, -mp.mpf(1)/2, -mp.mpf(1)/4, -mp.mpf(1)/5, mp.mpf(1)/16]]),
  '1/3': (mp.mpf(1)/3, [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log3^2', mp.log(3)**2), ('log23', mp.log(2)*mp.log(3))],
          [(str(a), a) for a in [mp.mpf(1)/2, mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/4, mp.mpf(3)/4, -mp.mpf(1)/2, -mp.mpf(1)/3, -mp.mpf(1)/4, mp.mpf(1)/9, mp.mpf(3)/8]]),
  '4/5': (mp.mpf(4)/5, [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log5^2', mp.log(5)**2), ('log25', mp.log(2)*mp.log(5))],
          [(str(a), a) for a in [mp.mpf(1)/2, mp.mpf(1)/5, mp.mpf(2)/5, mp.mpf(3)/5, mp.mpf(4)/5, mp.mpf(1)/9, mp.mpf(4)/9, -mp.mpf(1)/2, -mp.mpf(1)/5, -mp.mpf(1)/9, mp.mpf(5)/9]]),
}
for nm, (c, logs, li2s) in cases.items():
    T = T_alpha(c)
    names = ['T'] + [n for n, _ in logs] + ['Li2(' + n + ')' for n, _ in li2s] + ['1']
    v = [T] + [x for _, x in logs] + [mp.polylog(2, a) for _, a in li2s] + [mp.mpf(1)]
    rel = mp.pslq(v, tol=mp.mpf(10)**(-58), maxcoeff=500, maxsteps=1000)
    print(f"c={nm}:")
    if rel is None:
        print("   none")
        continue
    # normalize so coeff of T is -1
    iT = 0
    if rel[0] == 0:
        print("   T not in relation:", [(cf, names[i]) for i, cf in enumerate(rel) if cf]); continue
    scale = -rel[0]
    terms = [(mp.nstr(mp.mpf(cf)/scale, 12), names[i]) for i, cf in enumerate(rel) if cf != 0]
    print("   T = ", [(cf, nm2) for cf, nm2 in terms if nm2 != 'T'])
    s = sum((mp.mpf(cf)/scale)*v[i] for i, cf in enumerate(rel))
    print("   residual:", mp.nstr(abs(s), 3))
print("STEP13 DONE")
