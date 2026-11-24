"""M4f: verify an exported R5 form (planner/R5_nice_reduced.json format).
 1. EXACT full point symbol of classical part == target_symbol() (= symbol(W5)-symbol(N5)).
 2a. every Li5 atom (incl. Li5(b^4), Li5(1-b^4)) reduces to 0 in the exact (3,2) cobracket.
 2b. reduce(N5 - target(W5)) = 0 exactly.
 3. FRESH value at ~1700 digits (no cached values): R5 + N5 - W5_5400; expect ~1e-1596.
Precision set AFTER all imports.
"""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
from flint import arb, ctx
import mpmath as mp
import w5_symbol as W
import w5_cobracket as C
from w5_cobracket import reduce as cred, delta_word, target as ctarget
from w5_search import evaluate
import nice_search as NS

ZERO, ONE, Q, R, A, B = W.ZERO, W.ONE, W.Q, W.R, W.A, W.B
def log(s): print(s, flush=True)

EXTRA_Z = {'Li5((-2+2*sqrt2))': -2+2*R, 'Li5((17-12*sqrt2))': B**4,
           'Li5((1/2-1/4*sqrt2))': Q(F(1,2))-R/4, 'Li5((1/2+1/4*sqrt2))': Q(F(1,2))+R/4,
           'Li5((1/2-1/2*sqrt2))': Q(F(1,2))-R/2, 'Li5((1/2-3/8*sqrt2))': Q(F(1,2))-3*R/8,
           'Li5((-16+12*sqrt2))': -16+12*R}

def qnum(p,q):
    return mp.mpf(p.numerator)/p.denominator + mp.sqrt(2)*mp.mpf(q.numerator)/q.denominator

def main(path, digits=1620):
    terms = json.loads((HERE/path).read_text())['terms']
    sol = NS.solve(list(NS.S0), 'vjson')
    recmap = {nm: rec for nm,_,rec in sol['basis']}
    symmap = {nm: s for nm,s,_ in sol['basis']}
    for nm in terms:
        if nm in EXTRA_Z:
            z = EXTRA_Z[nm]; symmap[nm] = W.li_symbol(5,z); recmap[nm] = ('li',5,z)
        elif nm in ('z5','Lchi8_5'):
            symmap[nm] = {}; recmap[nm] = None
    co = {nm: (F(p),F(q)) for nm,(p,q) in terms.items()}
    # 1. EXACT SYMBOL
    t = {}
    for nm,(a,b) in co.items():
        if a or b:
            t = W.add(t, W.scale(symmap[nm], Q(a,b)))
    ts = NS.target_symbol()
    assert t == ts, f'{path}: SYMBOL MISMATCH'
    log(f'1. EXACT full point symbol == target_symbol: PASS')
    # 2. EXACT COBRACKET
    C.ban_float()
    for nm,(a,b) in co.items():
        if nm.startswith('Li5') and (a or b):
            z = recmap[nm][2]
            assert not cred(delta_word(ZERO,(ZERO,ZERO,ZERO,ZERO,ONE/z),ONE)).d, nm
    log(f'2a. every Li5 atom (incl. b^4, 1-b^4) reduces to 0 in (3,2) cobracket: PASS')
    cert = json.loads((HERE/'cobracket.json').read_text())
    right = C.P()
    for term in cert['correction']:
        x,y,c = (Q(*map(F,term[k])) for k in ['x','y','coefficient'])
        word = (ZERO,ZERO,ONE/x,ZERO,ONE/(x*y))
        right += delta_word(ZERO, word, ONE)*c
    delta = cred(right - ctarget())
    assert not delta.d
    log(f'2b. EXACT (3,2) cobracket: reduce(N5 - target(W5)) = 0: PASS')
    # 3. FRESH VALUE
    ctx.prec = math.ceil((digits+140)*math.log2(10)); mp.mp.dps = digits+120
    DD = digits+60
    total = mp.mpf(0)
    for nm,(a,b) in co.items():
        if not (a or b): continue
        if nm == 'z5':
            v = mp.mpf(arb(5).zeta().mid().str(DD+100, radius=False))
        elif nm == 'Lchi8_5':
            v = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
        else:
            v = mp.mpf(evaluate(recmap[nm], DD).mid().str(DD+100, radius=False))
        total += qnum(a,b)*v
    n5v = mp.mpf(0)
    for term in cert['correction']:
        x,y,c = (Q(*map(F,term[k])) for k in ['x','y','coefficient'])
        v = evaluate(('mpl',3,2,x,y), DD)
        n5v += qnum(c.p,c.q)*mp.mpf(v.mid().str(DD+100, radius=False))
    saved = json.loads((HERE/'W5_5400.json').read_text())
    res = total + n5v - mp.mpf(saved['value'])
    log(f'3. {path}: R5 + N5 - W5(5400 bits) = {mp.nstr(res, 8)}  (fresh {digits}-digit evaluation)')

if __name__=='__main__':
    for p in sys.argv[1:] or ['R5_best.json']:
        main(p)
