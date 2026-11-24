
import sys, json, re, pickle
sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr

# my S_F
d = json.load(open('../../10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
def parse(term):
    m = re.match(r'c\^(-?\d+)\s*(.*)', term[0]); e = int(m.group(1)); rest = m.group(2).strip()
    w = [int(x) for x in rest[2:-1].split(',')] if rest.startswith('H') else None
    return e, w
gmap = {0: C, 1: vneg(P), -1: Q}
S_F = {}
for t in d:
    e, w = parse(t)
    if w is None or len(w)!=4: continue
    coeff = Fr(t[1])
    k = tuple(gmap[a] for a in reversed(w))
    cur = S_F.setdefault(k, {})
    cur[e] = cur.get(e, Fr(0)) + coeff
    if cur[e]==0: del cur[e]

# COR4 F_symbol.json
fs = json.load(open('../COR4_M3alpha_classical/F_symbol.json'))
def parse_vec(s):
    # e.g. "00-10" -> parse as ints: sign-separated digits; digits are 0,1 or -1
    # format: string like "00-10" meaning (0,0,-1,0): parse char by char: '-' attaches to next digit
    out = []
    i = 0
    while i < len(s):
        if s[i] == '-':
            out.append(-int(s[i+1])); i += 2
        else:
            out.append(int(s[i])); i += 1
    return tuple(out)
S_F_cor4 = {}
for key, val in fs.items():
    m = re.match(r'c\^(-?\d+)\s+(.*)', key)
    e = int(m.group(1))
    vecs = tuple(parse_vec(v) for v in m.group(2).split('|'))
    cur = S_F_cor4.setdefault(vecs, {})
    cur[e] = cur.get(e, Fr(0)) + Fr(val)
    if cur[e]==0: del cur[e]

print("my tensors:", len(S_F), " cor4 tensors:", len(S_F_cor4))
# compare
keys = set(S_F) | set(S_F_cor4)
diff = 0
for k in keys:
    a = S_F.get(k, {}); b = S_F_cor4.get(k, {})
    if a != b:
        diff += 1
        if diff <= 12:
            print("DIFF at", k, " mine:", {e:str(x) for e,x in a.items()}, " cor4:", {e:str(x) for e,x in b.items()})
print("total differing tensors:", diff)
