script = r'''"""M4e: emit closed_form_nice.md for R5_nice (lchi_final_1200.json), grouped."""
import json
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
HERE = Path(__file__).resolve().parent

d = json.loads((HERE/'lchi_final_1200.json').read_text())
names = d['names']; ca=[F(c) for c in d['ca']]; cb=[F(c) for c in d['cb']]
aL=F(d['aL']); bL=F(d['bL'])
cert = json.loads((HERE/'cobracket.json').read_text())

def q2s(a, b):
    """latex of a + b*sqrt2 as a coefficient (with parens if needed)."""
    if b == 0:
        return sp.latex(sp.Rational(a.numerator, a.denominator)), False
    if a == 0:
        return sp.latex(sp.Rational(b.numerator, b.denominator)*sp.sqrt(2)), False
    s = f'{sp.latex(sp.Rational(a.numerator,a.denominator))}{"+" if b>0 else "-"}{sp.latex(sp.Rational(abs(b).numerator,abs(b).denominator))}\\sqrt2'
    return s, True

def arg_latex(nm):
    # extract argument inside Li5(...) / Li4(...) / Li3(...)
    z = nm[nm.index('(')+1:nm.rindex(')')]
    z = z.replace('*sqrt2', '\\sqrt2')
    return z

def atom_tex(nm):
    if nm.startswith('Li5('): return rf'\operatorname{{Li}}_5\!\left({arg_latex(nm)}\right)'
    if nm.startswith('Li4('): return rf'\operatorname{{Li}}_4\!\left({arg_latex(nm)}\right)'
    if nm.startswith('Li3('): return rf'\operatorname{{Li}}_3\!\left({arg_latex(nm)}\right)'
    # product: factors joined by *
    parts = []
    from collections import Counter
    for f in nm.split('*'):
        if f=='lam': parts.append('\\lambda')
        elif f=='L': parts.append('L')
        elif f=='pi2': parts.append('\\pi^2')
        elif f=='z3': parts.append('\\zeta(3)')
        elif f.startswith('Li2('): parts.append(rf'\operatorname{{Li}}_2\!\left({arg_latex(f)}\right)')
        elif f.startswith('Li3('): parts.append(rf'\operatorname{{Li}}_3\!\left({arg_latex(f)}\right)')
        elif f.startswith('Li4('): parts.append(rf'\operatorname{{Li}}_4\!\left({arg_latex(f)}\right)')
    # combine powers
    c = Counter(parts)
    out=[]
    for p in dict.fromkeys(parts):
        n=c[p]
        if n>1:
            if p in ('\\lambda','L'): out.append(p+rf'^{n}')
            elif p=='\\pi^2': out.append(rf'\pi^{{{2*n}}}')
            else: out.append(rf'\left[{p}\right]^{{{n}}}')
        else: out.append(p)
    return '\\,'.join(out)

def term_line(a, b, tex, count):
    s, paren = q2s(a, b)
    if paren: s = rf'\left({s}\right)'
    sign = '-' if (a<0 or (a==0 and b<0)) else ('+' if count else '')
    if sign=='-': s = q2s(-a, -b)[0];  s = rf'\left({s}\right)' if paren else s
    return sign + (s if s!='1' else '') + tex + r'\\ &'

groups = {'Li5':[], 'Li4':[], 'w32':[], 'pure':[]}
for nm,a,b in zip(names,ca,cb):
    if not (a or b): continue
    if nm.startswith('Li5('): groups['Li5'].append((nm,a,b))
    elif 'Li4(' in nm: groups['Li4'].append((nm,a,b))
    elif ('Li3' in nm) or ('z3' in nm) or ('Li2' in nm): groups['w32'].append((nm,a,b))
    else: groups['pure'].append((nm,a,b))

lines = ['# W5 = R5_nice: readable closed form with the weight-5 constant $L(\\chi_8,5)$','',
 'Let $b=\\sqrt2-1$, $\\lambda=\\log2$, $L=\\log(1+\\sqrt2)$, and',
 '$L(\\chi_8,5)=\\sum_{n\\ge1}\\chi_8(n)/n^5$ with $\\chi_8(n)=+1$ ($n\\equiv\\pm1\\bmod8$), $-1$ ($n\\equiv\\pm3\\bmod8$).',
 'The convention is $\\operatorname{Li}_{s,t}(x,y)=\\sum_{n>m>0}x^ny^m/(n^sm^t)$.','',
 'Numerical identity, verified to 1596 digits against the certified W5; its full point symbol and',
 '(3,2) Lie cobracket agree exactly with W5. Not an exact identity certificate.','']

def emit(title, items, start_count):
    lines.append(rf'\[ \text{{({title})}} \]') if False else None
    lines.append(r'\['); lines.append(r'\begin{aligned}')
    lines.append(rf'\text{{{title}}}\\')
    c = start_count
    for nm,a,b in items:
        lines.append(term_line(a,b,atom_tex(nm),c)); c+=1
    lines.append(r'\end{aligned}'); lines.append(r'\]')
    return c

# N5 first
lines.append(r'\['); lines.append(r'\begin{aligned}')
lines.append(r'\text{depth-two } \operatorname{Li}_{3,2}\text{ terms } (N_5, \text{ exact cobracket}):\\')
c=0
for term in cert['correction']:
    x,y,cf = (F(*map(str,term[k])) if k!='coefficient' else None for k in [])  # placeholder
for term in cert['correction']:
    x=Q_(term['x']); y=Q_(term['y']); cf=Q_(term['coefficient'])
    a,b = cf
    tex = rf'\operatorname{{Li}}_{{3,2}}\!\left({str(x).replace("*sqrt2",chr(92)+"sqrt2")},{str(y).replace("*sqrt2",chr(92)+"sqrt2")}\right)'
    lines.append(term_line(a,b,tex,c)); c+=1
lines.append(r'\end{aligned}'); lines.append(r'\]')
lines.append('')

lines.append(r'\['); lines.append(r'\begin{aligned}')
lines.append(r'L(\chi_8,5)\text{ and }\zeta(5)\text{ terms}\\')
s,paren=q2s(aL,bL)
lines.append((r'\left(' + s + r'\right)' if paren else s) + r'L(\chi_8,5)')
lines.append(r'\end{aligned}'); lines.append(r'\]')
lines.append('')

for title,key in [('Li_5 terms','Li5'),('\\operatorname{Li}_4\\cdot(\\text{weight 1})','Li4'),
                  ('(\\text{weight 3})\\cdot(\\text{weight 2})','w32'),('pure products','pure')]:
    emit(title, groups[key], 0)
    lines.append('')

(HERE/'closed_form_nice.md').write_text('\n'.join(lines)+'\n')
print('wrote closed_form_nice.md')
print(f'Li5 {len(groups["Li5"])}, Li4 {len(groups["Li4"])}, w32 {len(groups["w32"])}, pure {len(groups["pure"])}')
'''
# need Q_ helper for QS2 from string pair
helper = '''
import sys
sys.dont_write_bytecode=True
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
from qsqrt2 import QS2 as Q
from fractions import Fraction as F
def Q_(pair):
    return Q(F(pair[0]), F(pair[1]))
'''
M4E = '<project>/research/08_theory/M4e_w5_lchi'
open(M4E+'/display_nice.py','w').write(helper + script)
print("written display_nice.py")
