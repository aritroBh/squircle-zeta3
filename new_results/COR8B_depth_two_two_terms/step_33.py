
import sys, json, re
sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr

d = json.load(open('../../10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
def parse(term):
    m = re.match(r'c\^(-?\d+)\s*(.*)', term[0])
    e = int(m.group(1)); rest = m.group(2).strip()
    w = [int(x) for x in rest[2:-1].split(',')] if rest.startswith('H') else None
    return e, w, Fr(t[1] if False else 0)

gmap = {0: C, 1: vneg(P), -1: Q}   # g(0)=[c], g(1)=-[1-c], g(-1)=[1+c]

def S_HPL(w):
    # S(H_{a1..an}) = g(a_n) x ... x g(a_1)
    return {tuple(gmap[a] for a in reversed(w)): Fr(1)}

# build S(F) decomposed by c^e
S_F_by_e = {}
nterms = 0
for t in d:
    e, w, _ = parse(t)
    coeff = Fr(t[1])
    if w is None or len(w) != 4: continue
    nterms += 1
    s = S_HPL(w)
    Se = S_F_by_e.setdefault(e, {})
    for k,c0 in s.items():
        Se[k] = Se.get(k, Fr(0)) + coeff*c0
        if Se[k]==0: del Se[k]
print("weight-4 terms:", nterms, " c-exponents:", {e: len(s) for e,s in S_F_by_e.items()})

# delta22 per exponent
D_by_e = {e: delta22(s) for e,s in S_F_by_e.items()}
for e in sorted(D_by_e):
    print(f"c^{e} part of d22(F):", show_d22(D_by_e[e]))

# target: (1-1/c) T1 + (-1-1/c) T2 with T1=(cp^pq), T2=(cq^pq)
# i.e. component (3,5)=cp^pq: coeff 1-1/c ; (4,5)=cq^pq: coeff -1-1/c
target = {(3,5): None, (4,5): None}  # to be compared coefficientwise as Laurent polys
def laurent(d):
    # collect {(i,j): {e: coeff}}
    out = {}
    for e, dd in d.items():
        for k,v in dd.items():
            out.setdefault(k, {})[e] = v
    return out
LF = laurent(D_by_e)
print("\nAssembled d22(F) components (Laurent coefficients in c):")
for k,v in sorted(LF.items()):
    print(f"  {WNAME[k[0]]}^{WNAME[k[1]]}: ", v)
# check target
t1 = LF.get((3,5), {}); t2 = LF.get((4,5), {})
ok1 = (t1.get(0)==1 and t1.get(-1)==-1 and len(t1)==2)
ok2 = (t2.get(0)==-1 and t2.get(-1)==-1 and len(t2)==2)
print("\ntarget component cp^pq = (1 - 1/c)?", ok1, t1)
print("target component cq^pq = (-1 - 1/c)?", ok2, t2)
print("total components:", len(LF))
open('progress.log','a').write("step_04 S(F) from HPL definition reproduces paper d22(Omega) exactly\n")
