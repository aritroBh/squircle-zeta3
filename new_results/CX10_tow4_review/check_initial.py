"""Symbolic cache/Bell bridge and exact initial-tower samples."""
from pathlib import Path
import pickle
import sys
from fractions import Fraction as F
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'TOW4_tower_depth8'))
import sympy as sp
import tow2lib as T2
import tow3lib as T3
import tow4lib as T4
from check_data import Q, make_eval
d, m = T4.d, T4.m
for r in range(1, 9):
    assert sp.cancel(sp.together(T2.L(r)-T4._L[r])) == 0, r
    print('PASS L-cache == proved L(r), symbolic r =', r, flush=True)
for n in range(9):
    assert sp.simplify(sp.expand(T2.tau1(n)-T4.tau1f(n))) == 0, n
    assert sp.simplify(sp.expand(T2.tau2(n)-T4.tau2f(n))) == 0, n
    assert sp.conjugate(T4.tau1f(n)) == (-1)**n*T4.tau2f(n), n
print('PASS tau definitions and conjugation n=0..8; tau8 =', T4.tau1f(8), T4.tau2f(8), flush=True)
# Independent coefficient recurrence for exp(sum L_r y^r), rather than partitions.
a = [sp.Integer(1)]
for n in range(1, 9):
    a.append(sp.cancel(sum(r*T4._L[r]*a[n-r] for r in range(1, n+1))/n))
    assert sp.cancel(a[n]-T4.arho(n)) == 0, n
    assert sp.cancel(a[n]-T4.forced(n, 0)) == 0, n
print('PASS Bell recurrence == partition assembly == F(delta,0), delta=0..8', flush=True)
fresh = pickle.load(open(Path(__file__).resolve().parent/'forced_F78.pkl', 'rb'))
original = pickle.load(open(ROOT/'TOW4_tower_depth8'/'forced_F78.pkl', 'rb'))
assert fresh == original
print('PASS fresh and submitted forced caches identical, entries =', len(fresh), flush=True)
for dl in range(7):
    for s in range(dl//2+1):
        assert sp.cancel(T4.forced(dl, s)-T3.e(dl, s)) == 0, (dl, s)
for dl in range(9):
    assert sp.cancel(T4.forced(dl, dl//2)-T2.e_top(dl)) == 0, dl
print('PASS F == TOW3 e at all 16 indices delta<=6; F(top)==zigzag delta<=8', flush=True)
print('F(7,3) =', sp.factor(T4.forced(7, 3)), flush=True)
print('F(8,4) =', sp.factor(T4.forced(8, 4)), flush=True)
H = pickle.load(open(ROOT/'TOWN1_tower_large_n'/'data'/'towers_D8_K52.pkl', 'rb'))
for dl in range(9):
    ev = make_eval(T4.arho(dl), (d, m))
    count = diag = adjacent = 0
    for dd in (6, 7, 8, 12, 20, 40, 52):
        for pp in sorted({dd, dd-1, (dd+1)//2}):
            key = (pp, dd-pp)
            rho0 = Q.of(H[key][0][dd])
            assert not rho0.iszero()
            er, ei = ev(dd, pp)
            pred = rho0*Q(er, 0, ei, 0)
            assert (pred-Q.of(H[key][dl][dd])).iszero(), (key, dl)
            count += 1
            diag += int(key[0] == key[1])
            adjacent += int(key[0] == key[1]+1)
    print(f'PASS initial towers depth {dl}: {count} exact samples ({diag} diagonal, {adjacent} adjacent)', flush=True)
