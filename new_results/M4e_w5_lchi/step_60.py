clean = r'''"""M4e: emit closed_form_nice.md for R5_nice (lchi_final_1200.json), grouped."""
import json, sys
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
sys.dont_write_bytecode=True
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
from qsqrt2 import QS2 as Q
HERE = Path(__file__).resolve().parent

d = json.loads((HERE/'lchi_final_1200.json').read_text())
names = d['names']; ca=[F(c) for c in d['ca']]; cb=[F(c) for c in d['cb']]
aL=F(d['aL']); bL=F(d['bL'])
cert = json.loads((HERE/'cobracket.json').read_text())

def rat(x):
    return sp.latex(sp.Rational(x.numerator, x.denominator))
def q2s(a, b):
    if b == 0: return rat(a), False
    if a == 0: return sp.latex(sp.Rational(b.numerator,b.denominator)*sp.sqrt(2)), False
    s = f'{rat(a)}{"+" if b>0 else "-"}{rat(abs(b))}\\sqrt2'
    return s, True
def arg_latex(nm):
    z = nm[nm.index('(')+1:nm.rindex(')')].replace('*sqrt2','\\sqrt2')
    return z
def atom_tex(nm):
    for k in [5,4,3,2]:
        if nm.startswith(f'Li{k}('): return rf'\operatorname{{Li}}_{k}\!\left({arg_latex(nm)}\right)'
    from collections import Counter
    parts=[]
    for f in nm.split('*'):
        if f=='lam': parts.append('\\lambda')
        elif f=='L': parts.append('L')
        elif f=='pi2': parts.append('\\pi^2')
        elif f=='z3': parts.append('\\zeta(3)')
        elif f.startswith('Li2('): parts.append(rf'\operatorname{{Li}}_2\!\left({arg_latex(f)}\right)')
        elif f.startswith('Li3('): parts.append(rf'\operatorname{{Li}}_3\!\left({arg_latex(f)}\right)')
        elif f.startswith('Li4('): parts.append(rf'\operatorname{{Li}}_4\!\left({arg_latex(f)}\right)')
    c=Counter(parts); out=[]
    for p in dict.fromkeys(parts):
        n=c[p]
        if n>1:
            if p in ('\\lambda','L'): out.append(p+rf'^{n}')
            elif p=='\\pi^2': out.append(rf'\pi^{{{2*n}}}')
            else: out.append(rf'\left[{p}\right]^{{{n}}}')
        else: out.append(p)
    return '\\,'.join(out)
def term_line(a, b, tex, count):
    neg = (a<0) or (a==0 and b<0)
    aa, bb = (-a,-b) if neg else (a,b)
    s, paren = q2s(aa, bb)
    if paren: s = rf'\left({s}\right)'
    sign = '-' if neg else ('+' if count else '')
    return sign + ('' if s=='1' else s) + tex + r'\\ &'

groups = {'Li5':[], 'Li4':[], 'w32':[], 'pure':[]}
for nm,a,b in zip(names,ca,cb):
    if not (a or b): continue
    if nm.startswith('Li5('): groups['Li5'].append((nm,a,b))
    elif 'Li4(' in nm: groups['Li4'].append((nm,a,b))
    elif ('Li3' in nm) or ('z3' in nm) or ('Li2' in nm): groups['w32'].append((nm,a,b))
    else: groups['pure'].append((nm,a,b))

L=[]
L.append('# $W_5 = R_5^{\\rm nice}$: readable closed form with the weight-5 constant $L(\\chi_8,5)$\n')
L.append('Let $b=\\sqrt2-1$, $\\lambda=\\log2$, $L=\\log(1+\\sqrt2$, and '
         '$L(\\chi_8,5)=\\sum_{n\\ge1}\\chi_8(n)/n^5$ ($\\chi_8=+1$ on $n\\equiv\\pm1\\bmod8$, $-1$ on $\\pm3$). '
         '$\\operatorname{Li}_{s,t}(x,y)=\\sum_{n>m>0}x^ny^m/(n^sm^t)$.\n')
L.append('Numerical identity, verified to 1596 digits against the certified W5; its full point symbol and '
         '(3,2) Lie cobracket agree exactly with W5. Not an exact identity certificate.\n')
L.append('Total: ' + str(sum(len(groups[g]) for g in groups)+4+1) + ' terms '
         f'({len(groups["Li5"])} $\\operatorname{{Li}}_5$; 4 $\\operatorname{{Li}}_{{3,2}}$; 1 $L(\\chi_8,5)$; '
         f'{len(groups["Li4"])} $\\operatorname{{Li}}_4\\cdot$w1; {len(groups["w32"])} w3$\\cdot$w2; '
         f'{len(groups["pure"])} pure products).\n')

def block(title, items):
    L.append(r'\['); L.append(r'\begin{aligned}')
    L.append(rf'\text{{{title}}}\\')
    for i,(nm,a,b) in enumerate(items):
        L.append(term_line(a,b,atom_tex(nm),i))
    L[-1]=L[-1].removesuffix(r'\\ &')
    L.append(r'\end{aligned}'); L.append(r'\]'); L.append('')

# N5 block
L.append(r'\['); L.append(r'\begin{aligned}')
L.append(r'N_5\ (\operatorname{Li}_{3,2},\ \text{exact cobracket correction})\\')
for i,term in enumerate(cert['correction']):
    x=Q(F(term['x'][0]),F(term['x'][1])); y=Q(F(term['y'][0]),F(term['y'][1]))
    cf=Q(F(term['coefficient'][0]),F(term['coefficient'][1]))
    tex=rf'\operatorname{{Li}}_{{3,2}}\!\left({str(x).replace("*sqrt2",chr(92)+"sqrt2")},{str(y).replace("*sqrt2",chr(92)+"sqrt2")}\right)'
    L.append(term_line(cf.p,cf.q,tex,i))
L[-1]=L[-1].removesuffix(r'\\ &')
L.append(r'\end{aligned}'); L.append(r'\]'); L.append('')

# L(chi8,5) block
s,paren=q2s(aL,bL)
L.append(r'\[')
L.append(r'+\,\left(' + s + r'\right)L(\chi_8,5)' if paren else r'+\,'+s+r'L(\chi_8,5)')
L.append(r'\]'); L.append('')

block(r'\operatorname{Li}_5', groups['Li5'])
block(r'\operatorname{Li}_4\cdot(\text{weight }1)', groups['Li4'])
block(r'(\text{weight }3)\cdot(\text{weight }2)', groups['w32'])
block(r'\text{pure products}', groups['pure'])
(HERE/'closed_form_nice.md').write_text('\n'.join(L)+'\n')
print('wrote closed_form_nice.md, bytes:', len('\n'.join(L)))
'''
M4E = '<project>/research/08_theory/M4e_w5_lchi'
open(M4E+'/display_nice.py','w').write(clean)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice', PY, 'display_nice.py'], cwd=M4E, capture_output=True, text=True, timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
