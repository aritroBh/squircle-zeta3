"""r01: term-by-term exact comparison of the supplement's W5 table (paper/ieee/supplement_W5.tex, read only)
with M4f_w5_extras/R5_best.json (R_14) and M4e_w5_lchi/cobracket.json (N_5). Own parser; no project code."""
import re, json, sys
from fractions import Fraction as F
from collections import Counter
ROOT='<project>'
tex=open(ROOT+'/paper/ieee/supplement_W5.tex').read()
body=tex.split(r'\section*{The classical part}')[1].split(r'\section*{Status}')[0]
nbody=tex.split(r'\section*{The depth-two part $N_5$}')[1].split(r'\section*{The classical part}')[0]
def num(s):
    s=s.strip()
    m=re.fullmatch(r'\\tfrac\{(\d+)\}\{(\d+)\}',s)
    if m: return F(int(m.group(1)),int(m.group(2)))
    return F(int(s))
def coef(s):
    """returns (p,q) from a tex coefficient like '-\\bigl(561+\\tfrac{2695}{4}\\sqrt2\\bigr)' or '+1608' or '+\\tfrac{867}{16}\\sqrt2' """
    s=s.strip(); sg=1
    if s.startswith('-'): sg=-1; s=s[1:]
    elif s.startswith('+'): s=s[1:]
    if s.startswith(r'\bigl('):
        inner=s[len(r'\bigl('):-len(r'\bigr)')]
        m=re.fullmatch(r'(-?)(\\tfrac\{\d+\}\{\d+\}|\d+)([+-])(\\tfrac\{\d+\}\{\d+\}|\d+)\\sqrt2',inner)
        assert m, inner
        p=num(m.group(2))*(-1 if m.group(1) else 1); q=num(m.group(4))*(1 if m.group(3)=='+' else -1)
        return (sg*p,sg*q)
    if s.endswith(r'\sqrt2'): return (F(0),sg*num(s[:-len(r'\sqrt2')]))
    if s=='' : return (F(sg),F(0))
    return (sg*num(s),F(0))
ARG={'b':(-1,1),'-b':(1,-1),'-b^{2}':(-3,2),r'\tfrac{1}{\sqrt2}':(0,F(1,2)),r'\tfrac12':(F(1,2),0),
     r'1-\tfrac{1}{\sqrt2}':(1,F(-1,2)),r'2-\sqrt2':(2,-1),'2b':(-2,2),r'\tfrac{2-\sqrt2}{4}':(F(1,2),F(-1,4)),
     r'\tfrac{2+\sqrt2}{4}':(F(1,2),F(1,4)),r'-\tfrac{b}{2}':(F(1,2),F(-1,2)),r'\tfrac{4-3\sqrt2}{8}':(F(1,2),F(-3,8)),
     'b^{4}':(17,-12),'1-b^{4}':(-16,12)}
def parse_atoms(s):
    out=[]
    s=s.replace(r'\,',' ').strip()
    toks=re.findall(r'\\operatorname\{Li\}_\{?(\d)\}?\\bigl\((.+?)\\bigr\)(\^\{2\})?|\\lambda(\^\{(\d)\})?|(?<![a-zA-Z\\])L(\^\{(\d)\})?|\\pi(\^\{(\d)\})|\\zeta\((\d)\)|L\(\\chi_\{8\},5\)',s)
    for t in re.finditer(r'\\operatorname\{Li\}_\{?(\d)\}?\\bigl\((.+?)\\bigr\)(\^\{2\})?|\\lambda(\^\{(\d)\})?|L\(\\chi_\{8\},5\)|L(\^\{(\d)\})?|\\pi\^\{(\d)\}|\\zeta\((\d)\)',s):
        g=t.group(0)
        if g.startswith(r'\operatorname'):
            n=int(t.group(1)); x=ARG[t.group(2)]; k=2 if t.group(3) else 1
            out+= [('Li',n,(F(x[0]),F(x[1])))]*k
        elif g.startswith(r'\lambda'): out+=[('lam',)]*int(t.group(5) or 1)
        elif g.startswith(r'L(\chi'): out+=[('Lchi',)]
        elif g.startswith('L'): out+=[('L',)]*int(t.group(7) or 1)
        elif g.startswith(r'\pi'): out+=[('pi2',)]*(int(t.group(8))//2)
        elif g.startswith(r'\zeta'): out+=[('z3',) if t.group(9)=='3' else ('z5',)]
    return tuple(sorted(out))
lines=[l for l in body.splitlines() if l.startswith('&') or l.startswith('W_5-N_5={}&')]
tab={}
for l in lines:
    l=l.split('&',1)[1].rstrip('\\').strip()
    # split coefficient and atoms: coefficient is up to the first '\,' at top level
    if l.startswith(r'\bigl(') or l.startswith(r'-\bigl(') or l.startswith(r'+\bigl('):
        i=l.index(r'\bigr)')+len(r'\bigr)'); c=coef(l[:i]); rest=l[i:]
    else:
        i=l.index(r'\,'); c=coef(l[:i]); rest=l[i:]
    k=parse_atoms(rest)
    assert k and k not in tab, l
    tab[k]=c
print('table classical terms:',len(tab))
# JSON (own parse of key names)
def jarg(s):
    s=s.strip('()')
    if 'sqrt2' not in s: return (F(s),F(0))
    s=s[:-len('*sqrt2')]
    i=max(s.rfind('+'),s.rfind('-'))
    if i<=0: return (F(0),F(s))
    return (F(s[:i]),F(s[i:]))
def jkey(name):
    out=[]
    for a in re.split(r'\*(?![^()]*\))',name):
        if a in ('lam','L','pi2','z3','z5'): out.append((a,))
        elif a=='Lchi8_5': out.append(('Lchi',))
        elif a=='Li2(b)': out.append(('Li',2,(F(-1),F(1))))
        else:
            m=re.fullmatch(r'Li(\d)\((.*)\)',a); out.append(('Li',int(m.group(1)),jarg(m.group(2))))
    return tuple(sorted(out))
J={jkey(k):(F(p),F(q)) for k,(p,q) in json.load(open(ROOT+'/research/08_theory/M4f_w5_extras/R5_best.json'))['terms'].items()}
print('JSON classical terms:',len(J))
miss=set(tab)^set(J); print('keys only in one:',miss)
bad=[(k,tab[k],J[k]) for k in tab if k in J and tab[k]!=J[k]]
print('coefficient mismatches:',len(bad))
for b in bad: print('  ',b)
# N5
nt=[l for l in nbody.splitlines() if 'Li}_{3,2}' in l]
NJ=json.load(open(ROOT+'/research/08_theory/M4e_w5_lchi/cobracket.json'))['correction']
AR={'b':(F(-1),F(1)),r'\tfrac{1}{\sqrt2}':(F(0),F(1,2)),'-b':(F(1),F(-1))}
NT={}
for l in nt:
    l=l.split('&',1)[1].rstrip('\\').strip()
    m=re.match(r'(.*?)\\,?\\operatorname\{Li\}_\{3,2\}\\bigl\((.+),(-?1)\\bigr\)',l)
    c=coef(m.group(1)); NT[(AR[m.group(2)],(F(int(m.group(3))),F(0)))]=c
NJd={((F(t['x'][0]),F(t['x'][1])),(F(t['y'][0]),F(t['y'][1]))):(F(t['coefficient'][0]),F(t['coefficient'][1])) for t in NJ}
print('N5 terms table/json:',len(NT),len(NJd),' identical:',NT==NJd)
print('TOTAL table terms (classical + N5):',len(tab)+len(NT))
