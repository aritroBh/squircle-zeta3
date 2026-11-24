import os, sys, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, '.')
BASE = '<project>'
OUT = os.path.join(BASE, 'paper/squircle_note/outputs')
import mpmath as mp
mp.mp.dps = 1250
M3file = mp.mpf(open(os.path.join(OUT,'M3_methodA_1200.txt')).read().strip())
tails = pickle.load(open('m3_tails_exact.pkl','rb'))
W3p_red = tails['W3p'].num()
T = tails['T'].num()
W4_from_file = M3file - W3p_red - T
print('W4 (paper route, 1200-digit M3) =')
print(mp.nstr(W4_from_file, 1200))

# theorem value at 1250 dps
sys.path.insert(0, os.path.join(BASE,'research/08_theory/A3_W4_proof'))
from w4_symbol import conjecture_terms
r2 = mp.sqrt(2)
def qnum(c):
    return mp.mpf(c.p.numerator)/c.p.denominator + r2*mp.mpf(c.q.numerator)/c.q.denominator
def atomval(a):
    if a[0]=='li':
        v = mp.polylog(a[1], qnum(a[2]))
        assert abs(mp.im(v)) < mp.mpf(10)**(-1200), a
        return mp.re(v)
    if a[0]=='log': return mp.log(qnum(a[1]))
    if a[0]=='pi2': return mp.pi**2
    if a[0]=='z3': return mp.zeta(3)
    raise ValueError(a)
W4th = mp.mpf(0)
for c, atoms in conjecture_terms():
    t = qnum(c)
    for a in atoms: t *= atomval(a)
    W4th += t
print('W4 (theorem) =')
print(mp.nstr(W4th, 1200))
print('difference =', mp.nstr(W4_from_file - W4th, 8))
