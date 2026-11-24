"""Shared: load the two forms exactly from project JSON, parse atom names into Q(sqrt2) arguments."""
import json, re
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
HERE=Path(__file__).resolve().parent
TH=HERE.parent
R14_PATH=TH/'M4f_w5_extras/R5_best.json'     # paper's fourteen-argument form (70 classical terms; supplement table)
R5_PATH =TH/'M4f_w5_extras/R5_min3.json'     # the certified three-extra form (66 terms; W5C1/W5C2/W5G2)
N5_PATH =TH/'M4e_w5_lchi/cobracket.json'
W5_PATH =TH/'M4f_w5_extras/W5_5400.json'
def load(path):
    t=json.load(open(path))['terms']
    return {k:(F(p),F(q)) for k,(p,q) in t.items()}
def parse_arg(s):
    """'(1/2-3/8*sqrt2)' -> (Fraction p, Fraction q) exactly."""
    e=sp.nsimplify(sp.sympify(s.replace('sqrt2','sqrt(2)')))
    e=sp.expand(e); q=e.coeff(sp.sqrt(2)); p=sp.expand(e-q*sp.sqrt(2))
    assert p.is_rational and q.is_rational, s
    return F(int(p.p),int(p.q)), F(int(q.p),int(q.q))
def split_prod(name):
    parts=[]; depth=0; cur=''
    for ch in name:
        if ch=='(': depth+=1
        if ch==')': depth-=1
        if ch=='*' and depth==0: parts.append(cur); cur=''
        else: cur+=ch
    parts.append(cur); return parts
def atom_key(a):
    """atom string -> ('lam',)|('L',)|('pi2',)|('z3',)|('z5',)|('Lchi',)|('Li',n,(p,q))"""
    if a in ('lam','L','pi2','z3','z5'): return (a,)
    if a=='Lchi8_5': return ('Lchi',)
    if a=='Li2(b)': return ('Li',2,(F(-1),F(1)))
    m=re.fullmatch(r'Li(\d)\((.+)\)',a)
    assert m, a
    return ('Li',int(m.group(1)),parse_arg(m.group(2)))
def terms_keys(form):
    out={}
    for name,c in form.items():
        key=tuple(sorted(atom_key(a) for a in split_prod(name)))
        assert key not in out
        out[key]=c
    return out
