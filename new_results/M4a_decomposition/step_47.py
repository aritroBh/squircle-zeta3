import sys
sys.dont_write_bytecode = True
sys.path.insert(0,'.')
sys.path.insert(0, '<project>/research/08_theory/A3_W4_proof')
import mpmath as mp
from flint import acb, ctx
from w4_symbol import conjecture_terms
from qsqrt2 import QS2

mp.mp.dps = 1250
ctx.prec = 4200
r2 = mp.sqrt(2)
def qnum_mp(c): return mp.mpf(c.p.numerator)/c.p.denominator + r2*mp.mpf(c.q.numerator)/c.q.denominator
def qnum_fl(c): return acb(c.p.numerator)/c.p.denominator + acb(2).sqrt()*acb(c.q.numerator)/c.q.denominator
def atom_mp(a):
    if a[0]=='li': return mp.re(mp.polylog(a[1], qnum_mp(a[2])))
    if a[0]=='log': return mp.log(qnum_mp(a[1]))
    if a[0]=='pi2': return mp.pi**2
    if a[0]=='z3': return mp.zeta(3)
def atom_fl(a):
    if a[0]=='li': return acb(mp.nstr(qnum_mp(a[2]),1200)).polylog(acb(a[1]))
    if a[0]=='log': return qnum_fl(a[1]).log()
    if a[0]=='pi2': return acb.pi()**2
    if a[0]=='z3': return acb(3).zeta()

for idx,(c,atoms) in enumerate(conjecture_terms()):
    t_mp = qnum_mp(c)
    for a in atoms: t_mp *= atom_mp(a)
    t_fl = qnum_fl(c)
    for a in atoms: t_fl *= atom_fl(a)
    d = abs(t_mp - mp.mpf(t_fl.real.mid().str(1200,radius=False,more=True)))
    if d > mp.mpf('1e-1100'):
        print(f'term {idx}: coeff={c} atoms={atoms}  |mp-fl| = {mp.nstr(d,4)}')
print('per-term scan done')
