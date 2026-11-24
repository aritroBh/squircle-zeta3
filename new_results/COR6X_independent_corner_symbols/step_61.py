t = r'''
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
from cor6x_core import *
import sympy as sp
x = sp.Symbol('x')
print("S(G(1,0; x)) =", sym_G((1,0), x), "  (expect [x][1-x] = (x, 1-x))")
print("S(G(1,0; 1/2)) =", sym_expanded(sym_G((1,0), sp.Rational(1,2))), " (expect {(2,2): 1})")
print("S(G(0,1; 1/2)) =", sym_expanded(sym_G((0,1), sp.Rational(1,2))), " (expect {(2,2): -1})")
print("S(G(1,-1; 1/3)) =", sym_expanded(sym_G((1,-1), sp.Rational(1,3))))
print("  hand: [1+x]x[1-x] - [1+x]x[2] + [1-x]x[2] at x=1/3:")
print("   = [4/3][2/3] - [4/3][2] + [2/3][2] =", )
import itertools
# manual: (2*2-3)(2-3) etc.
def letters_Q(S):
    return S
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols")
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])
