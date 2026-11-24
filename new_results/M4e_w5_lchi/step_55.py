script = r'''"""M4e: full verification of R5_nice (lchi_final_{D}.json).

R5_nice = N5 + sum_n (ca[n]+sqrt2 cb[n]) atom(n) + (aL+sqrt2 bL) L(chi8,5).
Checks:
 1. EXACT symbol of (R5_nice - N5) equals target_symbol() (= symbol(W5) - symbol(N5)).
 2. EXACT (3,2) cobracket: classical part reduces to 0 (each Li5 -> 0), and
    reduce(N5_cobracket - target()) = 0  [M4c's proven identity, re-run].
 3. Fresh value of R5_nice minus W5_5400 at high precision.
 4. R5_nice minus M4c R5 (closed_form.json) at 300 digits.
"""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
import flint
from flint import arb, ctx
import mpmath as mp
import w5_symbol as W
import w5_cobracket as C
from w5_cobracket import reduce as cred, delta_word, target as ctarget, depth2
from w5_search import evaluate, mpl, numeric
import nice_search as NS

HERE = Path(__file__).resolve().parent
ZERO, ONE, Q, R, A, B = W.ZERO, W.ONE, W.Q, W.R, W.A, W.B
def log(s): print(s, flush=True)

def qnumf(p, q):
    return mp.mpf(p.numerator)/p.denominator + mp.sqrt(2)*mp.mpf(q.numerator)/q.denominator

def main(D=1200):
    d = json.loads((HERE/f'lchi_final_{D}.json').read_text())
    names = d['names']; ca = [F(c) for c in d['ca']]; cb = [F(c) for c in d['cb']]
    aL = F(d['aL']); bL = F(d['bL'])
    # recipes for the 70 atoms
    sol = NS.solve(list(NS.S0), 'verify')
    basis = sol['basis']; s0names = sol['names']
    assert names[:len(s0names)] == s0names
    symmap = {nm: s for nm, s, _ in basis}
    recmap = {nm: rec for nm, _, rec in basis}
    zmap = {'Li5((-2+2*sqrt2))': -2+2*R, 'Li5((17-12*sqrt2))': B**4,
            'Li5((1/2-1/4*sqrt2))': Q(F(1,2))-R/4, 'Li5((1/2+1/4*sqrt2))': Q(F(1,2))+R/4,
            'Li5((1/2-1/2*sqrt2))': Q(F(1,2))-R/2, 'Li5((1/2-3/8*sqrt2))': Q(F(1,2))-3*R/8}
    for nm in names[len(s0names):]:
        z = zmap[nm]
        symmap[nm] = W.li_symbol(5, z); recmap[nm] = ('li', 5, z)
    # 1. EXACT SYMBOL
    t = {}
    for nm, a, b in zip(names, ca, cb):
        if a or b:
            t = W.add(t, W.scale(symmap[nm], Q(a, b)))
    ts = NS.target_symbol()
    assert t == ts, 'symbol mismatch'
    log('1. EXACT full point symbol of R5_nice - N5 equals target_symbol (symbol(W5)-symbol(N5)): PASS')
    # 2. EXACT cobracket
    C.ban_float()
    for nm in names:
        if recmap[nm][0] == 'li' and recmap[nm][1] == 5 and (ca[names.index(nm)] or cb[names.index(nm)]):
            z = recmap[nm][2]
            assert not cred(delta_word(ZERO,(ZERO,ZERO,ZERO,ZERO,ONE/z),ONE)).d, nm
    log('2a. every Li5 atom in R5_nice reduces to 0 in the (3,2) cobracket: PASS')
    cert = json.loads((HERE/'cobracket.json').read_text())
    right = C.P()
    for term in cert['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x','y','coefficient'])
        word = (ZERO,ZERO,ONE/x,ZERO,ONE/(x*y))
        right += delta_word(ZERO, word, ONE)*c
    delta = cred(right - ctarget())
    assert not delta.d
    log('2b. EXACT (3,2) cobracket: reduce(N5 - target(W5)) = 0; and classical part = 0: PASS')
    # 3. VALUE vs W5_5400
    ctx.prec = math.ceil((1640)*math.log2(10)); mp.mp.dps = 1620
    total = mp.mpf(0)
    for nm, a, b in zip(names, ca, cb):
        if a or b:
            v = evaluate(recmap[nm], 1560)
            total += qnumf(a, b)*mp.mpf(v.mid().str(1600, radius=False))
    Lchi = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
    total += qnumf(aL, bL)*Lchi
    # N5
    n5v = mp.mpf(0)
    for term in cert['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x','y','coefficient'])
        v = evaluate(('mpl',3,2,x,y), 1560)
        n5v += qnumf(c.p, c.q)*mp.mpf(v.mid().str(1600, radius=False))
    total += n5v
    saved = json.loads((HERE/'W5_5400.json').read_text())
    res3 = total - mp.mpf(saved['value'])
    log(f'3. R5_nice - W5(5400 bits): {mp.nstr(res3, 10)}  (expect ~1e-1596)')
    # 4. vs M4c R5
    sys.path.insert(0, str(HERE.parent/'M4c_w5_closed_form'/'planner'))
    mp.mp.dps = 400; ctx.prec = math.ceil(460*math.log2(10))
    import planner_eval_R5 as E
    tot4 = mp.mpf(0)
    for nm, a, b in zip(names, ca, cb):
        if a or b:
            v = evaluate(recmap[nm], 340)
            tot4 += qnumf(a, b)*mp.mpf(v.mid().str(380, radius=False))
    tot4 += qnumf(aL, bL)*mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
    n5v4 = mp.mpf(0)
    for term in cert['correction']:
        x, y, c = (Q(*map(F, term[k])) for k in ['x','y','coefficient'])
        v = evaluate(('mpl',3,2,x,y), 340)
        n5v4 += qnumf(c.p, c.q)*mp.mpf(v.mid().str(380, radius=False))
    tot4 += n5v4
    res4 = tot4 - E.R5
    log(f'4. R5_nice - M4c R5: {mp.nstr(res4, 8)}  (expect < 1e-300)')
    log('ALL VERIFICATION CHECKS DONE')

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 1200)
'''
M4E = '<project>/research/08_theory/M4e_w5_lchi'
open(M4E+'/verify_nice.py','w').write(script)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice', PY, 'verify_nice.py', '1200'], cwd=M4E, capture_output=True, text=True, timeout=2700)
print(r.stdout); print("STDERR:", r.stderr[-2000:]); print(f"elapsed {time.monotonic()-t0:.0f}s")
open(M4E+'/progress.log','a').write("step_48: verify_nice.py 1200 (symbol, cobracket, values)\n")
