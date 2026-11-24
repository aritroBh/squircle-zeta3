"""WT5 step F: full point (3,2) Lie cobracket of E1, E2, E3 and exact depth-two corrections.
Adapted from M5b w5pp_cobracket.py. Controls: Li3/Li5 projection checks; W5'' target == cobracket5pp.json.
E3 control: reduced cobracket must be 0 (classical at symbol level)."""
import json, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import wt5_E_pairings as WS
import w5_symbol as W
import w5_cobracket as C
from w5_cobracket import reduce as cred, delta_word, depth2
from w4_right_coproduct import P, PE, ZERO, ONE, Q, R, A, B, ban_float

def log(s): print(s, flush=True)

def target_of(expr):
    ans = P()
    for (w, k), c in expr.items():
        if k == 0:
            ans += (delta_word(ONE, w, A) - delta_word(ONE, w, B)) * c
    return cred(ans)

CANDS = [(B, ONE), (R/2, ONE), (-B, ONE), (B, -ONE), (R/2, -ONE),
         (ONE-R/2, ONE), (2-R, ONE), (-B*B, ONE)]
ATOMS = [('z3',), (('Li3', B),), (('Li3', -B),), (('Li3', R/2),)]

def solve_depth2(lhs, label):
    import itertools
    reds = [(x, y, depth2(x, y)) for x, y in CANDS]
    for subset in itertools.combinations(range(len(reds)), 4):
        sel = [reds[i] for i in subset]
        if not all(c.q == 0 for _, _, red in sel for c in red.d.values()):
            continue
        matrix = flint.fmpq_mat([[str(red.d.get(at, ZERO).p) for _, _, red in sel] for at in ATOMS])
        if matrix.det() == 0:
            continue
        rhs = flint.fmpq_mat([[str(lhs.d.get(at, ZERO).p), str(lhs.d.get(at, ZERO).q)] for at in ATOMS])
        sol = matrix.inv() * rhs
        residual = lhs
        corr = []
        for i, (x, y, red) in enumerate(sel):
            c = Q(F(str(sol[i, 0])), F(str(sol[i, 1])))
            residual -= red * c
            corr.append((x, y, c))
        if not residual.d:
            log(f'EXACT solution for {label} with subset {subset}:')
            for x, y, c in corr:
                log(f'   {c} * Li32({x},{y})')
            return subset, corr
    return None, None

if __name__ == '__main__':
    ban_float()
    for z in [B, -B, R/2, ONE-R/2, 2-R]:
        red = cred(C.g3(ZERO, (ZERO, ZERO, ONE/z), ONE) + PE.LI3(z))
        assert not red.d, (z, red)
    log('Weight-three projection versus classical Li3 controls: PASS')
    assert not cred(delta_word(ZERO, (ZERO, ZERO, ZERO, ZERO, ONE/B), ONE)).d
    log('Classical Li5 cobracket control: PASS')
    # CONTROL: W5'' target and correction reproduce cobracket5pp.json
    t5 = target_of(WS.w5pp_expression())
    ref = json.loads((HERE/'cobracket5pp.json').read_text())
    log(f'W5pp target: {t5}')
    assert str(t5) == ref['target'], (str(t5), ref['target'])
    log('CONTROL: W5pp cobracket target == cobracket5pp.json: PASS')
    out = {'W5pp_target': str(t5)}
    for name, which in [('E1', 'D1'), ('E2', 'D2'), ('E3', 'D3')]:
        lhs = target_of(WS.e_expression(which))
        log(f'{name} full point delta32, tensor Li2(b): {lhs}')
        subset, corr = solve_depth2(lhs, name)
        assert corr is not None, f'no exact N({name}) found'
        out[name] = {'target': str(lhs), 'subset': list(subset),
                     'correction': [dict(x=[str(x.p), str(x.q)], y=[str(y.p), str(y.q)],
                                         coefficient=[str(c.p), str(c.q)]) for x, y, c in corr]}
    (HERE/'cobracket_e.json').write_text(json.dumps(out, indent=2) + '\n')
    with open(HERE/'progress.log', 'a') as f:
        f.write('step_F cobracket: Li3/Li5 controls PASS; W5pp target PASS; ' +
                '; '.join(f"{nm}: N={len(out[nm]['correction'])} atoms" for nm in ['E1','E2','E3']) + '\n')
    log('saved cobracket_e.json')
