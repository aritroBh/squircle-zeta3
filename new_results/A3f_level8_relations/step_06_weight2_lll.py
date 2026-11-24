"""step_06_weight2_lll.py: LLL probes at weight 2.
(a) all 88 known rows evaluate to ~0 at 330 digits (pipeline check)
(b) LLL on the 10 free-word values -> the missing relation (hint)
(c) probes: Li2(b) elementary? Cl2(pi/4) vs G?
"""
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
import mpmath as mp
mp.mp.dps = 335
from flint import fmpz_mat, fmpz
from fractions import Fraction as F

data = json.loads(open('weight2_values.json').read())
cols = [tuple(w) for w in data['cols']]
VAL = {tuple(eval(k)): mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) for k, v in data['values'].items()}
free = [tuple(w) for w in data['free_words']]

# (a) check known rows
import pickle
rows = pickle.load(open('audit_weight2.pickle','rb'))['rows']
worst = mp.mpf(0)
for row in rows:
    s = sum((mp.mpf(str(F(c))) * VAL[w] for w, c in row), mp.mpc(0))
    worst = max(worst, abs(s))
print('known rows worst residual at 330 digits:', mp.nstr(worst, 8), flush=True)

def lll_relations(vecs, scale_exp=320):
    """vecs: list of mp.mpc. Returns short integer vectors with tiny residual."""
    n = len(vecs)
    C = mp.power(10, scale_exp)
    M = fmpz_mat(n, n + 2)
    for i, v in enumerate(vecs):
        M[i, i] = 1
        M[i, n] = fmpz(int(mp.floor(C * mp.re(v) + mp.mpf('0.5'))))
        M[i, n+1] = fmpz(int(mp.floor(C * mp.im(v) + mp.mpf('0.5'))))
    B = M.lll()
    out = []
    for i in range(n):
        row = [int(B[i, j]) for j in range(n)]
        tail = [int(B[i, n]), int(B[i, n+1])]
        res = sum((mp.mpf(c) * v for c, v in zip(row, vecs)), mp.mpc(0))
        out.append((row, tail, res))
    return out

# (b) main LLL on the 10 free values
vecs = [VAL[w] for w in free]
out = lll_relations(vecs)
print('\nLLL on 10 free values (row, tail, residual):')
norms = []
for row, tail, res in out:
    nrm = max(abs(c) for c in row)
    norms.append((nrm, row, tail, res))
norms.sort(key=lambda t: t[0])
for nrm, row, tail, res in norms[:4]:
    print(' maxcoeff=', nrm, ' tail=', tail, ' residual=', mp.nstr(abs(res), 6))
    print('   coeffs=', row)

# (c) probes
sq2 = mp.sqrt(2); b = sq2 - 1; a = 1 + sq2
li2b = mp.polylog(2, b)
G = mp.catalan
cl = mp.im(mp.polylog(2, mp.e ** (mp.pi * mp.j / 4)))
named_re = [mp.pi**2, mp.log(2)**2, mp.log(a)**2, mp.log(2)*mp.log(a), li2b]
print('\nprobe real [pi^2, log2^2, loga^2, log2*loga, Li2(b)]:')
for row, tail, res in lll_relations([mp.mpc(x, 0) for x in named_re]):
    if abs(res) < mp.mpf('1e-300'):
        print('  RELATION:', row, mp.nstr(abs(res), 6)); break
else:
    print('  no relation found (smallest residuals:)',
          sorted(mp.nstr(abs(r),4) for _,_,r in [x for x in lll_relations([mp.mpc(x,0) for x in named_re])])[:3])
named_im = [G, cl, mp.pi*mp.log(2), mp.pi*mp.log(a)]
print('probe imag [G, Cl2(pi/4), pi*log2, pi*loga]:')
resid = lll_relations([mp.mpc(x, 0) for x in named_im])
for row, tail, res in resid:
    if abs(res) < mp.mpf('1e-300'):
        print('  RELATION:', row, mp.nstr(abs(res), 6)); break
else:
    print('  no relation; residuals:', sorted(mp.nstr(abs(r),4) for _,_,r in resid)[:3])
