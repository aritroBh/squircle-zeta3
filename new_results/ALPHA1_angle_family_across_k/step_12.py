"""step_12.py -- PSLQ identification of T(alpha) (weight two) at rational angles."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 70

def T_alpha(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

def basis_logs(args):
    out = {}
    for a in args:
        out[f'Li2({a})'] = mp.polylog(2, a)
    return out

cases = {
  '1/2': (mp.mpf(1)/2, [mp.pi**2, mp.log(2)**2, mp.log(3)**2, mp.log(2)*mp.log(3)],
          [mp.mpf(1)/2, mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/4, mp.mpf(3)/4, -mp.mpf(1)/2, -mp.mpf(1)/3, -mp.mpf(2)/3, -mp.mpf(1)/4, mp.mpf(1)/9]),
  '3/5': (mp.mpf(3)/5, [mp.pi**2, mp.log(2)**2, mp.log(5)**2, mp.log(2)*mp.log(5), mp.log(3)**2, mp.log(3)*mp.log(2), mp.log(3)*mp.log(5)],
          [mp.mpf(1)/2, mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/5, mp.mpf(2)/5, mp.mpf(3)/5, mp.mpf(4)/5, -mp.mpf(1)/2, -mp.mpf(1)/4, -mp.mpf(3)/4, -mp.mpf(1)/5, -mp.mpf(2)/5, -mp.mpf(3)/5, mp.mpf(1)/16]),
  '1/3': (mp.mpf(1)/3, [mp.pi**2, mp.log(2)**2, mp.log(3)**2, mp.log(2)*mp.log(3)],
          [mp.mpf(1)/2, mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/4, mp.mpf(3)/4, -mp.mpf(1)/2, -mp.mpf(1)/3, -mp.mpf(2)/3, -mp.mpf(1)/4, mp.mpf(1)/9, mp.mpf(1)/8, mp.mpf(3)/8]),
  '4/5': (mp.mpf(4)/5, [mp.pi**2, mp.log(2)**2, mp.log(5)**2, mp.log(2)*mp.log(5), mp.log(3)**2, mp.log(3)*mp.log(2), mp.log(3)*mp.log(5)],
          [mp.mpf(1)/2, mp.mpf(1)/5, mp.mpf(2)/5, mp.mpf(3)/5, mp.mpf(4)/5, mp.mpf(1)/9, mp.mpf(1)/3, mp.mpf(2)/3, -mp.mpf(1)/2, -mp.mpf(1)/5, -mp.mpf(2)/5, -mp.mpf(3)/5, -mp.mpf(1)/9, mp.mpf(4)/9, mp.mpf(5)/9]),
}
for nm, (c, logs, li2args) in cases.items():
    T = T_alpha(c)
    B = list(logs) + [mp.polylog(2, a) for a in li2args] + [mp.mpf(1)]
    names = (['pi^2'] + [f'log{x}' for x in logs[1:]] ) if False else None
    names = ['pi^2'] + [f'lb{i}' for i in range(len(logs)-1)] + [f'Li2({a})' for a in li2args] + ['1']
    v = [T] + B
    rel = mp.pslq(v, tol=mp.mpf(10)**(-60), maxcoeff=200, maxsteps=1000)
    print(f"c={nm}: T={mp.nstr(T,30)}")
    if rel is None:
        print("   no relation found")
    else:
        terms = [(coeff, names[i]) for i, coeff in enumerate(rel) if coeff != 0]
        print("   relation:", terms)
        # verify
        s = sum(coeff*v[i] for i, coeff in enumerate(rel))
        print("   residual:", mp.nstr(abs(s), 3))
print("STEP12 DONE")
