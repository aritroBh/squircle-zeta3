"""Re-use WT3's exact integration definitions for just the base and slope."""
import sys
from pathlib import Path
sys.dont_write_bytecode = True
source = Path(__file__).resolve().parent.parent / 'WT3_weight3_all_k' / 'step_02.py'
# The prefix defines integration/reduction rules and constants; the omitted suffix
# runs the k=3..15 batch and writes JSON outside this review directory.
prefix, marker, _ = source.read_text().partition('values = {}')
assert marker and 'with open(' not in prefix
exec(compile(prefix, str(source), 'exec'))

for name, R20, R11, expected in (
    ('slope', 1/(2*(v+1+S2))-1/(2*(v+S2-1)), 1/(2*v)-1/(v+1+S2), DELTA_rem),
    ('base', kernel_sympy(3,'20'), kernel_sympy(3,'11'), W3p_thm),
):
    print(f'Computing exact {name} integral', flush=True)
    raw = weight3_piece(2,0,R20)+weight3_piece(1,1,R11)
    clean, used, rank = PE.solve_li3_linear(PE.reduce_li2_all(raw), basis, extra_relations=KREL)
    assert not (clean-expected).d, name
    f20 = sp.lambdify(v,R20,modules='mpmath')
    f11 = sp.lambdify(v,R11,modules='mpmath')
    numerical = mp.quad(lambda x: f20(x)*lt(x)**2+f11(x)*lt(x)*l1t(x), [bN,1,aN])
    error = abs(expected.num()-numerical)
    assert error < mp.mpf('1e-60'), (name,error)
    print(f'PASS {name}: exact residual zero; 70-dps quadrature error {mp.nstr(error,6)}', flush=True)
