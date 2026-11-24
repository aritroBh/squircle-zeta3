"""Check printed formulas and Task D's finite algebra, including its missing premise."""
from pathlib import Path
import pickle
import re
import sys
sys.dont_write_bytecode = True
SOURCE = Path(__file__).resolve().parent.parent/'TOW4_tower_depth8'
sys.path.insert(0, str(SOURCE))
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor
import tow4lib as T4
d, m = T4.d, T4.m
text = (SOURCE/'proof_depth8.md').read_text()
local = {'d': d, 'm': m, 'i': sp.I, 'Pi_3': T4.Pi(3), 'Pi_4': T4.Pi(4)}
transforms = standard_transformations+(implicit_multiplication_application, convert_xor)
def parse(s):
    return parse_expr(s.replace('\n', ' '), local_dict=local, transformations=transforms)
errors = []
for dl, s, den in ((7, 1, 480), (7, 2, 6), (8, 1, 5760), (8, 2, 48), (8, 3, 4)):
    match = re.search(rf'  F\({dl},{s}\) = \((.*?)\) / \({den} Pi_{dl//2}\)', text, re.S)
    assert match, (dl, s)
    residual = sp.cancel(parse(match[1])/(den*T4.Pi(dl//2))-T4.forced(dl, s))
    if residual != 0:
        errors.append((dl, s))
    print('printed F', (dl, s), 'correct:', residual == 0, flush=True)
for dl in (7, 8):
    match = re.search(rf'  F\({dl},{dl//2}\) = (.*?)   \(zigzag', text)
    assert match
    residual = sp.cancel(parse(match[1])-T4.forced(dl, dl//2))
    if residual != 0:
        errors.append((dl, dl//2))
    print('printed top', dl, 'correct:', residual == 0, 'printed-minus-computed:', sp.factor(residual), flush=True)
assert errors == [(7, 3)], errors
for dl in (7, 8):
    saved = pickle.load(open(SOURCE/f'a{dl}_closed.pkl', 'rb'))
    expr = sp.sympify(saved, locals={'d': d, 'm': m})
    assert sp.cancel(expr-T4.arho(dl)) == 0
    print(f'PASS saved a{dl}_closed.pkl == Bell formula', flush=True)
match = re.search(r'N_7 =\s*(.*?)\n  a_8', text, re.S)
assert sp.cancel(parse(match[1])/(40320*T4.Pi(3))-T4.arho(7)) == 0
print('PASS printed a7 numerator == Bell formula', flush=True)

def feed(dl, values):
    return sum(sp.I*(2*d-1)/(2*d-dp-1)*(
        T4.tau1f(dl-dp)*values[dp].subs({d: d+1, m: m+1}, simultaneous=True)
        -T4.tau2f(dl-dp)*values[dp].subs({d: d+1, m: m}, simultaneous=True))
        for dp in range(dl+1))
for dl in (7, 8):
    P = dl//2
    U = {dp: sp.cancel(sum((-1)**s*T4.forced(dp, s) for s in range(dp//2+1))) for dp in range(dl+1)}
    master = sp.cancel(sp.together(feed(dl, U)-T4.arho(dl)))
    last = {dp: T4.forced(dp, P) for dp in range(dl+1)}
    closing = sp.cancel(feed(dl, last)-T4.forced(dl, P))
    assert master == closing == 0
    assert sp.cancel(closing-(-1)**P*master) == 0
    print(f'PASS Task D finite alternating identity and closing at depth {dl}', flush=True)

# Top equations alone do not fix the prescribed zigzag numerator. At depth 8,
# perturb only a8 by d^4/Pi4. Its forced perturbations are Delta_d^s(d^4)/Pi4.
poly = d**4
perturb = []
for s in range(6):
    perturb.append(sp.cancel(poly/T4.Pi(4)))
    poly = sp.expand(poly.subs(d, d+1)-poly)
for s in range(5):
    H = (2*d-1)/(2*d-9)*perturb[s].subs(d, d+1)
    assert sp.cancel(H-perturb[s]-perturb[s+1]) == 0
assert perturb[5] == 0 and sp.cancel(perturb[4]-24/T4.Pi(4)) == 0
assert perturb[1] != 0
print('PASS missing-premise example: a8 += d^4/Pi4 preserves closing and lower depths, but changes top by 24/Pi4', flush=True)
print('Thus closing alone cannot identify a separately prescribed zigzag-descending family; top matching must be supplied.', flush=True)
