"""corner_equals_squircle.py -- NUMERICAL check of the (hand-proved) identity M(u; pi/2) = -F(u)/2 (squircle).
Proof (result.md Sec. 1): substitute t = g/(1+g) in F = 2 int_0^{1/2} Phi dt; the difference of the two integrands
integrates to (1 - 2^-u)/u - (1 - 2^-u)/u = 0, using int_0^1 g^(u-1) (1+g)^(-1-u) dg = 2^-u/u.
Also counts the weight-4 HPL words of B3 by prefactor."""
import sys
import mpmath as mp
sys.path.insert(0, '.')
from check_b3_split import R3
mp.mp.dps = 40
for u in (mp.mpf('0.3'), mp.mpf('0.7'), mp.mpf('-0.4')):
    corner = mp.quad(lambda g: ((1 + g)**(-1 - u)*mp.sqrt(1 + g**(2 - 2*u)) - 1)*g**(u - 1), [0, mp.mpf(1)/8, 1])
    sq = -mp.quad(lambda t: t**(u - 1) + (1 - t)**(u - 1) - mp.sqrt(t**(2*u - 2) + (1 - t)**(2*u - 2)), [0, mp.mpf(1)/8, mp.mpf(1)/2, mp.mpf(7)/8, 1])/2
    print(f'u = {u}: corner M(u;pi/2) = {mp.nstr(corner, 25)}, squircle -F/2 = {mp.nstr(sq, 25)}, diff = {mp.nstr(corner - sq, 3)}')
ws = {}
for e, w, v in R3:
    if len(w) == 4:
        ws.setdefault(e, set()).add(w)
print('weight-4 HPL words in B3 by prefactor c^e:', {e: len(s) for e, s in sorted(ws.items())}, '; distinct words overall:', len(set().union(*ws.values())))
