"""a09_V_needed2.py -- as a08 but with strict admissibility (constant of 1-f is +-2^k) and a numerical
niceness test: value v at c=1/sqrt2 is 'nice' if v (or 1/v when v<0, by exact inversion) is one of
 1/sqrt2, -1/sqrt2, 1/2, 1-1/sqrt2, 2-sqrt2, b, -b, -b^2, b^2, b^4, -1, 0   (b = sqrt2-1),
i.e. reducible to the conjecture basis by D4/I4/V4 (Kimi_B's 11 reducible values); v > 1 allowed when 1/v is nice
(real part of Li4, exact inversion), provided the argument stays > 1 along the path (checked later)."""
import sys
exec(open('a07_V_feasibility.py').read().split("pt = 1/sp.sqrt(2)")[0])
import mpmath as mp
mp.mp.dps = 30
r2 = mp.sqrt(2); bb = r2 - 1
NV = [1/r2, -1/r2, mp.mpf(1)/2, 1 - 1/r2, 2 - r2, bb, -bb, -bb**2, bb**2, bb**4, mp.mpf(-1), mp.mpf(0)]
def isnice(v):
    if any(abs(v - x) < mp.mpf(10)**-20 for x in NV): return True
    if v < 0 and any(abs(1/v - x) < mp.mpf(10)**-20 for x in NV): return True
    if v > 1 and any(abs(1/v - x) < mp.mpf(10)**-20 for x in NV): return True   # Re Li4 + exact inversion
    return False
pt = 1/sp.sqrt(2)
cols = []
for f in adm:
    fv = sp.lambdify(c, f, 'mpmath')(1/r2)
    rs = rho(li4_symbol(f))
    if rs:
        cols.append((f, fv, rs))
target = rho(SQ)
keys = sorted(set(target) | set().union(*[set(r) for _, _, r in cols]))
M = sp.Matrix([[r.get(k, 0) for _, _, r in cols] for k in keys])
bvec = sp.Matrix([target.get(k, 0) for k in keys])
print('strictly admissible arguments:', len(adm), ' with nonzero rho:', len(cols))
print('rank all:', M.rank(), ' rank all+target:', M.row_join(bvec).rank())
nice = [i for i, (f, v, r) in enumerate(cols) if isnice(v)]
Mn = M[:, nice]
print('nice columns:', len(nice), ' rank', Mn.rank(), ' rank+target', Mn.row_join(bvec).rank())
non = [i for i in range(len(cols)) if i not in nice]
vals = {}
for i in non:
    vals.setdefault(mp.nstr(cols[i][1], 12), []).append(str(cols[i][0]))
print('non-nice values present (value: #args):')
for v, fs in sorted(vals.items(), key=lambda t: float(t[0])):
    print('   %16s : %d   e.g. %s' % (v, len(fs), fs[0]))
need = [i for i in non if M[:, nice + [i]].row_join(bvec).rank() == M[:, nice + [i]].rank()]
print('single non-nice columns completing the nice set:', [(str(cols[i][0]), mp.nstr(cols[i][1], 10)) for i in need])
