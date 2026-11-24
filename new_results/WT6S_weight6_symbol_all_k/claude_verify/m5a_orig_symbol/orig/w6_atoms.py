"""Survey of depth-two atom pairing vectors at weight 6 (for the N6 solve)."""
import json
import sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'A3_W4_proof'))
from w4_symbol import ban_float
import w6_symbol as W
from w4_symbol import Q, R, A, B, ONE as ONEQ, ZERO as ZEROQ

ban_float()

XS = {
    'b': B, '-b': -B, 'b^2': B*B, '-b^2': -B*B,
    '1/sqrt2': Q(F(0), F(1, 2)), '-1/sqrt2': Q(F(0), -F(1, 2)),
    '1/2': Q(F(1, 2)), '-1/2': Q(F(-1, 2)),
    '2': Q(2), '-2': Q(-2),
    '1-1/sqrt2': ONEQ - Q(F(0), F(1, 2)),
    '2-sqrt2': Q(2) - R,
}
YS = {'1': ONEQ, '-1': -ONEQ, '1/2': Q(F(1, 2)), '-1/2': Q(F(-1, 2)),
      '1/sqrt2': Q(F(0), F(1, 2)), 'b': B, '-b': -B, '2': Q(2), '-2': Q(-2)}

out = {}
for mn in [(4, 2), (3, 3)]:
    for xn, x in XS.items():
        for yn, y in YS.items():
            try:
                s = W.depth2_symbol(mn[0], mn[1], x, y)
                v = [W.pair33(s)] + W.pair42(s)
            except Exception as e:
                continue
            if all(t == ZEROQ for t in v):
                continue
            out[f"Li_{mn[0]}{mn[1]}({xn},{yn})"] = [str(t) for t in v]
print(f"atoms with nonzero pairing vector: {len(out)}")
for k, v in out.items():
    print(f"  {k}: (33, 42a, 42b, 42c) = {v}")
(HERE/'atom_pairings.json').write_text(json.dumps(out, indent=1) + '\n')
S6 = W.symbol(W.w6_expression())
print("W6: (33, 42a,42b,42c) =", str(W.pair33(S6)), [str(t) for t in W.pair42(S6)])
