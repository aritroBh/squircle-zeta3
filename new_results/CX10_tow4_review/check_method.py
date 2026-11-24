"""Compare slow/fast feed expressions in QQ(i)(p), including a nonzero control."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'TOW4_tower_depth8'))
import sympy as sp
import prove_adj56 as slow
import prove_adj78_fast as fast

field = sp.QQ_I.frac_field(fast.p)
for dl, s in ((5, 1), (5, 2), (6, 2), (6, 3)):
    lhs_s = slow.eR(dl, s, slow.d, slow.p)+slow.eR(dl, s+1, slow.d, slow.p)
    lhs_f = fast.eR(dl, s, fast.d, fast.p)+fast.eR(dl, s+1, fast.d, fast.p)
    assert field.from_sympy(lhs_s) == field.from_sympy(lhs_f)
    ts, tf = slow.feed_terms(dl, s), fast.feed_terms(dl, s)
    assert len(ts) == len(tf) == 3*(dl+1)
    assert all(field.from_sympy(a) == field.from_sympy(b) for a, b in zip(ts, tf))
    residual = field.from_sympy(lhs_f)-sum(map(field.from_sympy, tf), field.zero)
    assert residual == field.zero
    assert residual+field.one != field.zero
    print(f'PASS ({dl},{s}): lhs and every feed term slow==fast; field residual zero; +1 control nonzero', flush=True)
for dl in (7, 8):
    # Normalized E3 coefficients, not the adjacent child's normalization.
    import prove_E3b78 as e3
    top = sp.factor(e3.ehat(dl, dl//2))
    print(f'stable diagonal top depth {dl}, /rho_pp(0): {top}', flush=True)
print('PASS all slow/fast expression comparisons', flush=True)
