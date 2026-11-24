
import sys, time, pickle, itertools
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 170
import rul1, rul1num

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

l2 = mp.log(2)
z = {n: mp.zeta(n) for n in range(2, 11)}
Li = {n: mp.polylog(n, mp.mpf(1)/2) for n in range(4, 11)}
val = pickle.load(open('basis_vals_89.pkl','rb'))

prim_vals = [('log2', 1, l2), ('z2', 2, z[2]), ('z3', 3, z[3]), ('Li4', 4, Li[4]), ('z5', 5, z[5]), ('Li5', 5, Li[5]),
             ('Li6', 6, Li[6]), ('X6', 6, val['X6']), ('z7', 7, z[7]), ('Li7', 7, Li[7]), ('X7a', 7, val['X7a']),
             ('X7b', 7, val['X7b']), ('Li8', 8, Li[8]), ('X8a', 8, val['X8a']), ('X8b', 8, val['X8b']),
             ('X8c', 8, val['X8c']), ('X8d', 8, val['X8d'])]
def monomial_vals(W, maxw):
    prims = [(nm, wt, v) for nm, wt, v in prim_vals if wt <= maxw]
    out = []
    def rec(rem, start, acc, prod):
        if rem == 0:
            out.append(("*".join(prims[i][0] for i in acc), prod)); return
        for idx in range(start, len(prims)):
            nm, wt, v = prims[idx]
            if wt > rem: continue
            acc.append(idx); rec(rem-wt, idx, acc, prod*v); acc.pop()
    rec(W, 0, [], mp.mpf(1))
    return out

def E1(w, dps=160):
    return rul1num.g_word_over_0m1_num(w, dps=dps)
def I_val(m, j):
    return mp.factorial(j)*mp.factorial(m)*E1((0,)*(m+1)+(-1,)*j)

basis8 = monomial_vals(8, 7) + [('Li8', Li[8]), ('X8a', val['X8a']), ('X8b', val['X8b']), ('X8c', val['X8c']), ('X8d', val['X8d'])]
print("basis8:", len(basis8))
basis9 = monomial_vals(9, 8) + [('z9', z[9]), ('Li9', Li[9]), ('zbar(7,1,1)', val['E9_711']), ('zbar(6,2,1)', val['E9_621']),
                                ('zbar(5,3,1)', val['E9_531']), ('zbar(5,1,1,1,1)', val['E9_51111']), ('zbar(4,2,1,1,1)', val['E9_42111']),
                                ('X9fp', val['X9fp'])]
print("basis9:", len(basis9))

reds89 = {}
for n, bvals in [(8, basis8), (9, basis9)]:
    for m in range(n-1):
        j = n-1-m
        Iv = I_val(m, j)
        vec = [Iv] + [v for _, v in bvals]
        rel = mp.pslq(vec, tol=mp.mpf('1e-140'), maxcoeff=10**15, maxsteps=10000)
        if rel is None:
            print(f"I_{{{m},{j}}} (w{n}): PSLQ FAIL", flush=True)
            continue
        a0 = rel[0]
        coeffs = [Fr(-x, a0) for x in rel[1:]]
        res = Iv - sum(c*v for c, (_, v) in zip(coeffs, bvals))
        reds89[(m,j)] = coeffs
        nterms = sum(1 for c in coeffs if c)
        print(f"I_{{{m},{j}}} (w{n}): {nterms} terms, |resid| {mp.nstr(abs(res),3)}", flush=True)
pickle.dump({'reds': reds89, 'basis8': [nm for nm,_ in basis8], 'basis9': [nm for nm,_ in basis9]}, open('I_reductions_w89.pkl','wb'))
log("PSLQ reductions for w8/w9 I's computed")
