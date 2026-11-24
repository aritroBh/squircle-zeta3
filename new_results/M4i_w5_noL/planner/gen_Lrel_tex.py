"""Planner: LaTeX display of M4i's L(chi8,5) relation, built only from L_relation.json with M4g's formatting helpers.
Every coefficient is q*sqrt2 with q rational, so sqrt2 is factored out. Writes Lchi8_relation.tex (standalone)."""
import sys
from fractions import Fraction
from pathlib import Path
GEN = Path(__file__).resolve().parents[2] / 'M4g_w5_tex' / 'gen_supplement.py'
src = GEN.read_text()
src = src[:src.index('def fifty(v):')]
sys.argv = ['gen', str(Path(__file__).resolve().parents[1] / 'L_relation.json'), '/dev/null']
import os; os.chdir(GEN.parent)
exec(compile(src, str(GEN), 'exec'))
assert all(Fraction(p) == 0 for p, q in terms.values())
def line(key):
    s, body = coeff_latex(Fraction(terms[key][1]), Fraction(0))
    return '&%s%s\\,%s' % (s, body, atom_latex(key))
order = LI5_ORDER + LI4_ORDER + W3W2_ORDER + PURE_ORDER
rows = [line(k) for k in order]
rows[0] = rows[0].replace('&+', '&', 1)
body = ('\\begin{align*}\n\\frac{L(\\chi_8,5)}{\\sqrt2}={}' + '\\\\\n'.join(rows) + '\n\\end{align*}\n')
doc = ('\\documentclass[11pt]{article}\n\\usepackage{amsmath,amssymb}\n\\usepackage[margin=1in]{geometry}\n'
       '\\allowdisplaybreaks\n\\begin{document}\n\\noindent With $b=\\sqrt2-1$, $\\lambda=\\log2$, $L=\\log(1+\\sqrt2)$ '
       'and $L(\\chi_8,5)=\\sum_{n\\ge1}\\chi_8(n)n^{-5}$, numerically (agreement to 1740 digits; the symbol and the '
       '$(3,2)$ cobracket of the right-hand side vanish exactly):\n' + body + '\\end{document}\n')
out = Path(__file__).with_name('Lchi8_relation.tex'); out.write_text(doc)
print('wrote', out.name, 'terms', len(order))
