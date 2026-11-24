WORK='<project>/research/08_theory/M4f_w5_extras'
src = r'''"""M4f: generate w5_closed_form.tex (standalone article) and m4_terms.json (all coefficient
data the tex is generated from; the M4 numeric check parses THE SAME file).

W5 = R5_best.json (70 terms) + N5 (cobracket.json, 4 Li_{3,2} terms).
W4'' = M4a w4pp_closed_form.json (proved; = M4b2 R''_nice).
W3^(4) = W3' - zeta(3)/8 - pi^2 lam/24 - lam^3/12 (W3' from M4c planner_assemble_M4.py).
T = (pi^2 + 6 lam^2)/24.
"""
import json, re, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
T = HERE.parent

# ---------- exact coefficient data ----------
r5 = json.loads((HERE/'R5_best.json').read_text())['terms']
cert = json.loads((HERE/'cobracket.json').read_text())
N5 = []
for term in cert['correction']:
    x = F(term['x'][0]), F(term['x'][1]); y = F(term['y'][0]), F(term['y'][1])
    c = F(term['coefficient'][0]), F(term['coefficient'][1])
    N5.append(dict(x=[str(x[0]),str(x[1])], y=[str(y[0]),str(y[1])], c=[str(c[0]),str(c[1])]))
w4 = json.loads((T/'M4a_decomposition'/'w4pp_closed_form.json').read_text())['closed_form_terms']
# parse W4'' into {atomkey: [p,q]}
def w4parse(terms):
    out = {}
    for cs, nm in terms:
        c = F(cs); parts = nm.split('*')
        q = F(0)
        if parts[0] == 'sqrt2': q, parts = c, parts[1:]; c = F(0)
        key = '*'.join(parts) if parts else '1'
        a,b = out.get(key, [F(0),F(0)])
        out[key] = [a+c, b+q]
    return out
W4 = w4parse(w4)
# W3^(4): exact data (rational part p, sqrt2 part q)
W34 = {
 'z3': [F(-5,16), F(7,32)],
 'Li3(b)': [F(1), F(1)],
 'Li3(-b)': [F(1), F(-1)],
 'Li3(1/sqrt2)': [F(0), F(-2)],
 'Li2(b)*L': [F(0), F(2)],
 'Li2(b)*log2': [F(0), F(-1)],
 'pi^2*L': [F(-1,8), F(1,12)],
 'pi^2*log2': [F(-1,12), F(-1,24)],
 'log2^3': [F(-1,24), F(1,24)],
 'log2^2*L': [F(0), F(1,4)],
 'log2*L^2': [F(-1,2), F(-1,2)],
 'L^3': [F(5,6), F(-1,6)],
}
Tt = {'pi^2': [F(1,24), F(0)], 'log2^2': [F(1,4), F(0)]}

json.dump(dict(R5=r5, N5=N5, W4pp={k:[str(a),str(b)] for k,(a,b) in W4.items()},
               W34={k:[str(a),str(b)] for k,(a,b) in W34.items()},
               T={k:[str(a),str(b)] for k,(a,b) in Tt.items()}),
          open(HERE/'m4_terms.json','w'), indent=1)

# ---------- LaTeX helpers ----------
def rat(x):
    x = F(x)
    if x.denominator == 1: return str(x.numerator)
    return rf'\tfrac{{{x.numerator}}}{{{x.denominator}}}'
def coeff_tex(a, b):
    """coefficient (a + b sqrt2) as latex; returns (tex, negative, paren)."""
    a, b = F(a), F(b)
    neg = (a < 0) or (a == 0 and b < 0)
    if neg: a, b = -a, -b
    if b == 0:
        s = rat(a); paren = False
    elif a == 0:
        s = (rat(b) if b != 1 else '') + r'\sqrt2'; paren = False
    else:
        s = f'{rat(a)}{"+" if b>0 else "-"}' + ('' if abs(b)==1 else rat(abs(b))) + r'\sqrt2'
        paren = True
    if neg and b != 0 and a != 0: pass
    return s, neg, paren

FIELD_RE = re.compile(r'^\(?([+-]?[\d/]+)?([+-])?([\d/]+)?\*sqrt2\)?$|^[\d/]+$|^-1$')
def parse_qs2(s):
    """Parse QS2 repr like '1/2*sqrt2', '(1/2-3/8*sqrt2)', '-1', '(17-12*sqrt2)'."""
    s = s.strip()
    if s == 'b': return r'b'
    t = s.strip('()')
    if '*sqrt2' in t:
        head, q = t.split('*sqrt2')
        q = q or '1'
        if head in ('', '+'): p, sign = F(0), 1
        else:
            m = re.match(r'^([+-]?[\d/]+)([+-])$', head)
            p, sign = F(m.group(1)), (1 if m.group(2)=='+' else -1)
        qf = sign*F(q)
        if p == 0: return (rf'{rat(qf)}\sqrt2' if qf != 1 else r'\sqrt2')
        return f'{rat(p)}{"+" if qf>0 else "-"}' + ('' if abs(qf)==1 else rat(abs(qf))) + r'\sqrt2'
    return rat(F(t))

def atom_tex(a):
    if a == 'log2' or a == 'lam': return r'\lambda'
    if a == 'L': return 'L'
    if a == 'pi^2' or a == 'pi2': return r'\pi^2'
    if a == 'z3': return r'\zeta(3)'
    if a == 'z5': return r'\zeta(5)'
    if a == 'Lchi8_5': return r'L(\chi_8,5)'
    if a == 'Li2(b)': return r'\operatorname{Li}_2\!\left(b\right)'
    if a == 'Li3(b)': return r'\operatorname{Li}_3\!\left(b\right)'
    if a == 'Li3(-b)': return r'\operatorname{Li}_3\!\left(-b\right)'
    if a == 'Li3(1/sqrt2)': return r'\operatorname{Li}_3\!\left(\tfrac{1}{\sqrt2}\right)'
    m = re.fullmatch(r'Li(\d+)\((.*)\)', a)
    if m:
        return rf'\operatorname{{Li}}_{m.group(1)}\!\left({parse_qs2(m.group(2))}\right)'
    m = re.fullmatch(r'(log2|L)\^(\d)', a)
    if m:
        base = r'\lambda' if m.group(1)=='log2' else 'L'
        return base + '^{' + m.group(2) + '}'
    raise ValueError(a)

def split_atoms(nm):
    out, depth, cur = [], 0, ''
    for ch in nm:
        if ch == '(' : depth += 1
        if ch == ')' : depth -= 1
        if ch == '*' and depth == 0: out.append(cur); cur = ''
        else: cur += ch
    out.append(cur); return out

def prod_tex(nm):
    from collections import Counter
    parts = split_atoms(nm)
    cnt = Counter(parts)
    out = []
    for a in parts:
        if cnt[a] == 0: continue
        n = cnt[a]; cnt[a] = 0
        t = atom_tex(a)
        if n > 1:
            if a in ('lam','log2','L'): t += f'^{n}'
            elif a in ('pi2','pi^2'): t = rf'\pi^{{{2*n}}}'
            else: t = rf'\left[{t}\right]^{{{n}}}'
        out.append(t)
    return r'\,'.join(out)

def block(title, items):
    """items: list of (a, b, tex). One aligned display."""
    L = [r'\[', r'\begin{aligned}']
    L.append(rf'\text{{{title}}}\\')
    for i,(a,b,tex) in enumerate(items):
        s, neg, paren = coeff_tex(a,b)
        if paren: s = rf'\left({s}\right)'
        sign = '-' if neg else ('+' if i else '')
        if s == '1' and tex != '1': s = ''
        L.append(f'&{sign}{s}{tex}' + (r'\\' if i < len(items)-1 else ''))
    L.append(r'\end{aligned}'); L.append(r'\]')
    return L

# ---------- group the R5 terms ----------
groups = {'Li5':[], 'Li4':[], 'w32':[], 'pure':[]}
z5c = lc = None
for nm,(p,q) in r5.items():
    if nm == 'z5': z5c = (p,q); continue
    if nm == 'Lchi8_5': lc = (p,q); continue
    if nm.startswith('Li5('): groups['Li5'].append((nm,p,q))
    elif 'Li4' in nm: groups['Li4'].append((nm,p,q))
    elif ('Li3' in nm) or ('z3' in nm) or ('Li2' in nm): groups['w32'].append((nm,p,q))
    else: groups['pure'].append((nm,p,q))

doc = []
doc.append(r'\documentclass[11pt]{article}')
doc.append(r'\usepackage{amsmath,amssymb}')
doc.append(r'\usepackage[margin=1in]{geometry}')
doc.append(r'\newcommand{\Li}{\operatorname{Li}}')
doc.append(r'\begin{document}')
doc.append(r'\section*{The moment $M_4$ in closed form}')
doc.append(r'Write $b=\sqrt2-1$, $\lambda=\log 2$, $L=\log(1+\sqrt2)$, '
           r'$\Li_{s,t}(x,y)=\sum_{n>m>0}x^ny^m/(n^sm^t)$, and '
           r'$L(\chi_8,5)=\sum_{n\ge1}\chi_8(n)/n^5$ with $\chi_8(n)=+1$ for $n\equiv\pm1\ (8)$, $-1$ for $n\equiv\pm3\ (8)$. '
           r'Then $M_4 = W_5 + W_4'' + W_3^{(4)} + T$, where $T=(\pi^2+6\lambda^2)/24$ and the weight-$5$, '
           r'weight-$4$ and weight-$3$ pieces are as follows. '
           r'The weight-$5$ identity $W_5=R_5$ below is verified to $1596$ digits against the certified value of $W_5$, '
           r'and its full point symbol and $(3,2)$ Lie cobracket agree with $W_5$ exactly; '
           r'the assembled $M_4$ agrees with the high-precision quadrature value to $10^{-800}$.')
doc.append('')
doc.append(r'\subsection*{The weight-5 piece $W_5 = R_5$}')
doc.append(rf'$W_5 = N_5 + (\text{{classical terms below}})$, with {len(r5)+4} terms in total: '
           rf'4 $\Li_{{3,2}}$ (in $N_5$), {len(groups["Li5"])} $\Li_5$, $\zeta(5)$, $L(\chi_8,5)$, '
           rf'{len(groups["Li4"])} $\Li_4\cdot(\text{{weight }}1)$, {len(groups["w32"])} $(\text{{weight }}3)\cdot(\text{{weight }}2)$, '
           rf'{len(groups["pure"])} pure products.')
# N5 block
items = []
for t in N5:
    x = parse_qs2(t['x'][0] + ('*sqrt2' if F(t['x'][1]) else '')) if True else ''
    def q2(pp,qq):
        p, q = F(pp), F(qq)
        if q == 0: return rat(p)
        if p == 0: return ('' if abs(q)==1 else rat(abs(q)))+r'\sqrt2' if q>0 else '-'+(('' if abs(q)==1 else rat(abs(q)))+r'\sqrt2')
        return f'{rat(p)}{"+" if q>0 else "-"}' + ('' if abs(q)==1 else rat(abs(q))) + r'\sqrt2'
    xt = q2(t['x'][0], t['x'][1]); yt = q2(t['y'][0], t['y'][1])
    items.append((t['c'][0], t['c'][1], rf'\Li_{{3,2}}\!\left({xt},{yt}\right)'))
doc += block(r'$N_5$ (exact $(3,2)$-cobracket correction)', items)
# Li5 + zeta5 + Lchi
items = [(p,q, atom_tex(nm)) for nm,p,q in groups['Li5']]
items.append((z5c[0], z5c[1], r'\zeta(5)'))
items.append((lc[0], lc[1], r'L(\chi_8,5)'))
doc += block(r'$\Li_5$, $\zeta(5)$, $L(\chi_8,5)$', items)
doc += block(r'$\Li_4\cdot(\lambda,\,L)$', [(p,q,prod_tex(nm)) for nm,p,q in groups['Li4']])
doc += block(r'$(\text{weight }3)\cdot(\text{weight }2)$', [(p,q,prod_tex(nm)) for nm,p,q in groups['w32']])
doc += block(r'pure products $(\lambda,\,L,\,\pi)$', [(p,q,prod_tex(nm)) for nm,p,q in groups['pure']])
# W4''
doc.append(r'\subsection*{The weight-4 piece $W_4''$ (proved)}')
g4 = {'Li4':[], 'w31':[], 'w22':[], 'pure':[]}
for k,(a,b) in W4.items():
    if 'Li4' in k: g4['Li4'].append((a,b,k))
    elif ('Li3' in k) or ('z3' in k): g4['w31'].append((a,b,k))
    elif 'Li2' in k: g4['w22'].append((a,b,k))
    else: g4['pure'].append((a,b,k))
order4 = ['log2*log2*log2*log2','log2*log2*log2*L','log2*log2*L*L','log2*L*L*L','L*L*L*L',
          'pi^2*log2*log2','pi^2*log2*L','pi^2*L*L','Li2(b)*log2*log2','Li2(b)*log2*L','Li2(b)*L*L',
          'pi^2*pi^2','pi^2*Li2(b)','Li2(b)*Li2(b)','z3*log2','z3*L','Li3(b)*log2','Li3(b)*L',
          'Li3(-b)*log2','Li3(-b)*L','Li3(1/sqrt2)*log2','Li3(1/sqrt2)*L']
items = []
for nm,p,q in g4['Li4']: items.append((p,q,prod_tex(nm)))
doc += block(r'$W_4'':\ \Li_4$ terms', items)
items = []
for k in order4:
    if k in W4:
        a,b = W4[k]; items.append((a,b,prod_tex(k)))
doc += block(r'$W_4''$: products', items)
# W34
doc.append(r'\subsection*{The weight-3 piece $W_3^{(4)}$ (proved)}')
doc.append(r'$W_3^{(4)} = W_3'' - \zeta(3)/8 - \pi^2\lambda/24 - \lambda^3/12$; with like terms combined:')
order3 = ['z3','Li3(b)','Li3(-b)','Li3(1/sqrt2)','Li2(b)*L','Li2(b)*log2','pi^2*L','pi^2*log2',
          'log2^3','log2^2*L','log2*L^2','L^3']
items = [(W34[k][0], W34[k][1], prod_tex(k)) for k in order3]
doc += block(r'$W_3^{(4)}$', items)
# T
doc.append(r'\subsection*{The weight-2 piece $T$ (proved)}')
doc.append(r'\[ T = \tfrac{1}{24}\,\pi^2 + \tfrac{1}{4}\,\lambda^2 = \frac{\pi^2+6\lambda^2}{24}. \]')
doc.append(r'\subsection*{Assembly}')
doc.append(r'\[ M_4 = W_5 + W_4'' + W_3^{(4)} + T, \]')
doc.append(r'with the four pieces displayed above. No like terms remain between different weights; '
           r'within $W_3^{(4)}$ and $W_4''$ like terms have been combined as shown.')
doc.append(r'\end{document}')
(HERE/'w5_closed_form.tex').write_text('\n'.join(doc)+'\n')
print('wrote w5_closed_form.tex:', len('\n'.join(doc)), 'bytes; terms:', len(r5)+4, '+ W4', len(W4), '+ W34', len(W34), '+ T', len(Tt))
'''
open(f'{WORK}/gen_tex.py','w').write(src)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','gen_tex.py'],
                   cwd=WORK, capture_output=True, text=True, timeout=300)
print(r.stdout, r.stderr[-1200:])
