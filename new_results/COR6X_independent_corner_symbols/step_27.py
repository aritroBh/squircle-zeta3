t = r'''
# step_05a.py -- numerical verification of the differential formula for I(a0;a1..an;a_{n+1})
import mpmath as mp
mp.mp.dps = 50

def I_num(params):
    a0, aN = params[0], params[-1]
    letters = params[1:-1]
    n = len(letters)
    def rec(i, lo):
        # integral of letters[i:] over [lo, ...] with outermost letter letters[i]
        if i == n - 1:
            f = lambda t: 1/(t - letters[i])
            return mp.quad(f, [lo, aN])
        def f(t, i=i):
            return rec(i + 1, lo) if False else None
        # iterate: I = int_{a0}^{aN} dt1/(t1-a1) * I_inner(letters[1:], up to t1)
        def integrand(t):
            return rec_inner(i + 1, a0, t) / (t - letters[i])
        return mp.quad(integrand, [a0, aN])
    def rec_inner(i, lo, hi):
        # integral of letters[i:] with upper limit hi, lower lo, outermost letter letters[i]
        if i == n - 1:
            return mp.quad(lambda t: 1/(t - letters[i]), [lo, hi])
        return mp.quad(lambda t: rec_inner(i + 1, lo, t) / (t - letters[i]), [lo, hi])
    return rec_inner(0, a0, aN)

# test at distinct generic parameters
import random
base = [mp.mpf(0), mp.mpf('0.31'), mp.mpf('-0.17'), mp.mpf('0.13'), mp.mpf('0.8')]
n = 3
val = I_num(base)
print("I =", mp.nstr(val, 20))

# formula: dI = sum_i I(omit a_i) [ dlog(prev_i - a_i) - dlog(next_i - a_i) ], cyclic:
# prev_1 = a_{n+1}, next_n = a_0
def Iomit(params, i):
    return I_num([params[0]] + [p for j, p in enumerate(params[1:-1]) if j != i - 1] + [params[-1]])

h = mp.mpf('1e-12')
ok = True
for k in range(1, n + 2):   # vary params[k] (skip a0 = index 0 for a moment)
    dp = [p for p in base]
    dp[k] += h
    dnum = (I_num(dp) - val) / h
    # formula prediction: sum_i Iomit * d/dparams[k] log((prev_i - a_i)/(next_i - a_i))
    pred = mp.mpf(0)
    for i in range(1, n + 1):
        ai = base[i]
        prev = base[n + 1] if i == 1 else base[i - 1]
        nxt = base[0] if i == n else base[i + 1]
        # log derivatives wrt param index k
        d = mp.mpf(0)
        if k == i or (i == 1 and k == n + 1) or (i == n and k == 0):
            pass
        # d/dp log(prev - ai): p=prev: +1/(prev-ai); p=ai: -1/(prev-ai)
        for idx, sign in ((k, None),):
            pass
        # compute explicitly
        term = mp.mpf(0)
        # from log(prev - ai)
        if k == i:
            term += -1/(prev - ai) - 1/(ai - nxt)
        if (i == 1 and k == n + 1) or (i > 1 and k == i - 1):
            term += 1/(prev - ai)
        if (i == n and k == 0) or (i < n and k == i + 1):
            term += -1/(ai - nxt)
        pred += Iomit(base, i) * term
    err = abs(dnum - pred)
    print(f"d/d a[{k}]: numeric {mp.nstr(dnum,12)}  formula {mp.nstr(pred,12)}  err {mp.nstr(err,3)}")
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_05a.py","w").write(t)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_05a.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols")
print(r.stdout); print("ERR:", (r.stderr or "")[-1500:])
