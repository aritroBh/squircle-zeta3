"""Full point (3,2) Lie cobracket of W5'' and the exact depth-two correction N5''.
Adapted from M4c w5_cobracket.py: same machinery, W5'' expression as target.
"""
import json, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import w5pp_symbol as WS
import w5_symbol as W
import w5_cobracket as C
from w5_cobracket import reduce as cred, delta_word, depth2
from w4_right_coproduct import P, PE, ZERO, ONE, Q, R, A, B, ban_float

def log(s): print(s, flush=True)

def target5pp():
    ans = P()
    for (w, k), c in WS.w5pp_expression().items():
        if k == 0:
            ans += (delta_word(ONE, w, A) - delta_word(ONE, w, B)) * c
    return cred(ans)

if __name__ == '__main__':
    ban_float()
    # controls from M4c (its own checks run on import of w5_cobracket? no: rerun key ones)
    for z in [B, -B, R/2, ONE-R/2, 2-R]:
        red = cred(C.g3(ZERO, (ZERO, ZERO, ONE/z), ONE) + PE.LI3(z))
        assert not red.d, (z, red)
    log('Weight-three projection versus classical Li3 controls: PASS')
    assert not cred(delta_word(ZERO, (ZERO, ZERO, ZERO, ZERO, ONE/B), ONE)).d
    log('Classical Li5 cobracket control: PASS')
    lhs = target5pp()
    log(f'Full point delta32(W5pp), tensor Li2(b): {lhs}')
    # candidate depth-two atoms; solve exact square system as M4c did
    cands = [(B, ONE), (R/2, ONE), (-B, ONE), (B, -ONE), (R/2, -ONE),
             (ONE-R/2, ONE), (2-R, ONE), (-B*B, ONE)]
    reds = [(x, y, depth2(x, y)) for x, y in cands]
    atoms = [('z3',), (('Li3', B),), (('Li3', -B),), (('Li3', R/2),)]
    # try to find a 4-subset solving exactly; report
    import itertools
    sol_found = None
    for subset in itertools.combinations(range(len(reds)), 4):
        sel = [reds[i] for i in subset]
        if not all(c.q == 0 for _, _, red in sel for c in red.d.values()):
            continue
        matrix = flint.fmpq_mat([[str(red.d.get(at, ZERO).p) for _, _, red in sel] for at in atoms])
        if matrix.det() == 0:
            continue
        rhs = flint.fmpq_mat([[str(lhs.d.get(at, ZERO).p), str(lhs.d.get(at, ZERO).q)] for at in atoms])
        sol = matrix.inv() * rhs
        residual = lhs
        corr = []
        for i, (x, y, red) in enumerate(sel):
            c = Q(F(str(sol[i, 0])), F(str(sol[i, 1])))
            residual -= red * c
            corr.append((x, y, c))
        if not residual.d:
            sol_found = (subset, corr)
            log(f'EXACT solution with subset {subset}:')
            for x, y, c in corr:
                log(f'   {c} * Li32({x},{y})')
            break
    assert sol_found is not None, 'no exact N5pp correction found among candidates'
    subset, corr = sol_found
    out = {'target': str(lhs),
           'correction': [dict(x=[str(x.p), str(x.q)], y=[str(y.p), str(y.q)],
                               coefficient=[str(c.p), str(c.q)]) for x, y, c in corr],
           'subset': list(subset)}
    (HERE/'cobracket5pp.json').write_text(json.dumps(out, indent=2) + '\n')
    log('saved cobracket5pp.json')
