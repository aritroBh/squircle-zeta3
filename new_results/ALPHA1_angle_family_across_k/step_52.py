step16 = r'''"""step_16.py -- identify T(r) with a fixed basis in r (or c), coefficients rational, stable across angles."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 60

def T_alpha_c(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

def basis(c):
    r = (1-c)/(1+c)
    L = lambda x: mp.log(x)
    A = lambda x: mp.polylog(2, x)
    out = [('pi^2', mp.pi**2)]
    out += [('Li2(c)', A(c)), ('Li2(-c)', A(-c))]
    out += [('Li2(r)', A(r)), ('Li2(-r)', A(-r))]
    out += [('Li2((1-c)/2)', A((1-c)/2)), ('Li2((1+c)/2)', A((1+c)/2))]
    out += [('Li2((1-r)/2)', A((1-r)/2)), ('Li2((1+r)/2)', A((1+r)/2))]
    out += [('Li2(c/(1+c))', A(c/(1+c))), ('Li2(-c/(1+c))', A(-c/(1+c)))]
    out += [('Li2(1-c)', A(1-c)), ('Li2(1-r)', A(1-r))]
    out += [('log2^2', L(2)**2), ('logc^2', L(c)**2), ('log1pc^2', L(1+c)**2), ('log1mc^2', L(1-c)**2)]
    out += [('logr^2', L(r)**2), ('log1pr^2', L(1+r)**2), ('log1mr^2', L(1-r)**2)]
    out += [('log2logc', L(2)*L(c)), ('log2log1pc', L(2)*L(1+c)), ('logclog1pc', L(c)*L(1+c))]
    out += [('log2logr', L(2)*L(r)), ('logrlog1pr', L(r)*L(1+r)), ('logrlog1mr', L(r)*L(1-r))]
    return out

angles = [mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, mp.mpf(2)/7, mp.mpf(5)/6, mp.mpf(1)/5]
Ts = [T_alpha_c(c) for c in angles]
B0 = basis(angles[0])
names = [n for n, _ in B0]
print("basis size:", len(names))

# PSLQ at first angle with iterative removal
def identify(T, B):
    B = list(B)
    while True:
        v = [T] + [x for _, x in B]
        rel = mp.pslq(v, tol=mp.mpf(10)**(-48), maxcoeff=2000, maxsteps=3000)
        if rel is None: return None
        if rel[0] != 0:
            return rel, [n for n, _ in B]
        idx = max(i for i in range(1, len(rel)) if rel[i] != 0)
        del B[idx-1]

out = identify(Ts[0], B0)
if out:
    rel, kept = out
    coef = [mp.mpf(x)/(-rel[0]) for x in rel]
    print("relation at c=1/2:")
    for i, cf in enumerate(coef):
        if cf != 0 and i >= 1:
            print(f"   {mp.nstr(cf,8)}  {kept[i-1]}")
    # verify at all angles
    for c, T in zip(angles, Ts):
        B = basis(c)
        bmap = {n: x for n, x in B}
        pred = sum(cf*bmap[kept[i-1]] for i, cf in enumerate(coef) if i >= 1 and cf != 0)
        print(f"   c={mp.nstr(c,8)}: T={mp.nstr(T,22)} pred={mp.nstr(pred,22)} diff={mp.nstr(abs(T-pred),3)}")
print("STEP16 DONE")
'''
open('step_16.py','w').write(step16)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_16.py'], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
