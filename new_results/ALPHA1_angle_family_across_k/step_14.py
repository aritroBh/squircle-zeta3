"""step_14.py -- PSLQ T(alpha) with iteratively independent basis."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 80

def T_alpha(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

def identify(T, B):
    """B: list of (name, value). Find T = combo."""
    B = list(B)
    while True:
        names = ['T'] + [n for n, _ in B]
        v = [T] + [x for _, x in B]
        rel = mp.pslq(v, tol=mp.mpf(10)**(-70), maxcoeff=1000, maxsteps=2000)
        if rel is None:
            return None
        if rel[0] != 0:
            scale = -rel[0]
            terms = [(rel[i], B[i-1][0]) for i in range(1, len(rel)) if rel[i] != 0]
            s = sum(rel[i]*v[i] for i in range(len(rel)))
            # T = sum rel[i]/scale * B[i-1]
            val = sum(mp.mpf(rel[i])/scale*B[i-1][1] for i in range(1, len(rel)) if rel[i]!=0)
            return [(str(mp.mpf(rel[i])/scale), B[i-1][0]) for i in range(1, len(rel)) if rel[i] != 0], abs(T - val)
        # drop last nonzero basis element
        idx = max(i for i in range(1, len(rel)) if rel[i] != 0)
        del B[idx-1]

def li2set(args):
    return [(f'Li2({mp.nstr(a,5)})', mp.polylog(2, a)) for a in args]

cases = {}
H = mp.mpf(1)/2
def L(*args):
    return list(args)
cases['1/2'] = (mp.mpf(1)/2,
    [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log3^2', mp.log(3)**2), ('log2log3', mp.log(2)*mp.log(3))],
    [H, mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/4, mp.mpf(3)/4, -H, -mp.mpf(1)/3, -mp.mpf(2)/3, -mp.mpf(1)/4, -mp.mpf(3)/4, mp.mpf(1)/9])
cases['3/5'] = (mp.mpf(3)/5,
    [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log5^2', mp.log(5)**2), ('log2log5', mp.log(2)*mp.log(5)), ('log3^2', mp.log(3)**2)],
    [H, mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/5, mp.mpf(2)/5, mp.mpf(3)/5, mp.mpf(4)/5, -H, -mp.mpf(1)/4, -mp.mpf(1)/5, -mp.mpf(2)/5, -mp.mpf(3)/5, -mp.mpf(4)/5, mp.mpf(1)/16, mp.mpf(3)/16, mp.mpf(5)/16, mp.mpf(9)/16])
cases['1/3'] = (mp.mpf(1)/3,
    [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log3^2', mp.log(3)**2), ('log2log3', mp.log(2)*mp.log(3))],
    [H, mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/4, mp.mpf(3)/4, -H, -mp.mpf(1)/3, -mp.mpf(2)/3, -mp.mpf(1)/4, -mp.mpf(3)/4, mp.mpf(1)/9, mp.mpf(1)/8, mp.mpf(3)/8, mp.mpf(1)/16])
cases['4/5'] = (mp.mpf(4)/5,
    [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log5^2', mp.log(5)**2), ('log2log5', mp.log(2)*mp.log(5)), ('log3^2', mp.log(3)**2)],
    [H, mp.mpf(1)/5, mp.mpf(2)/5, mp.mpf(3)/5, mp.mpf(4)/5, mp.mpf(1)/9, mp.mpf(4)/9, mp.mpf(5)/9, -H, -mp.mpf(1)/5, -mp.mpf(2)/5, -mp.mpf(3)/5, -mp.mpf(1)/9, -mp.mpf(4)/9, mp.mpf(1)/16, mp.mpf(1)/25, mp.mpf(1)/8, mp.mpf(3)/8])
for nm, (c, logs, li2args) in cases.items():
    T = T_alpha(c)
    B = logs + li2set(li2args)
    out = identify(T, B)
    print(f"c={nm}: T={mp.nstr(T,25)}")
    if out:
        terms, err = out
        print("   T =", terms, " resid:", mp.nstr(err, 3))
    else:
        print("   no identification")
print("STEP14 DONE")
