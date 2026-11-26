"""r05 (referee): planted errors against the claim's own all-n checker (s03_prove.check_part, run on copies).
Each mutant must make at least one of: even identity, odd identity, low orders  fail."""
import sys, os, io, contextlib, copy, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "claim_copy")); os.chdir(os.path.join(HERE, "claim_copy"))
sys.setrecursionlimit(20000)
import s03_prove as S
import g8_allN as G
from g8_zsum import lin

P5 = G.mu_side(5)
Z = lambda **kw: tuple(kw.get(f"z{w}", 0) for w in range(2, 6))
mon = {'1': Z(), 'z2': Z(z2=1), 'z2^2': Z(z2=2), 'z3': Z(z3=1), 'z4': Z(z4=1), 'z5': Z(z5=1), 'z2*z3': Z(z2=1, z3=1)}
raw = {nm: json.load(open(f"hpl_fit_k5_{nm}.json")) for nm in mon}

def run(part, lines, label):
    terms = S.parse(lines)
    ev = G.hpl_side(terms, 0); diff = lin((1, ev), (-1, P5[mon[part]]))
    odd = G.hpl_side(terms, 1)
    lo = G.low_order(terms, lambda N: Fr(0), lo=-6, hi=1)
    caught = bool(diff) or bool(odd) or bool(lo)
    print(f"{label:70s} even-fails={bool(diff)!s:5} odd-fails={bool(odd)!s:5} low-fails={bool(lo)!s:5} -> {'CAUGHT' if caught else 'NOT CAUGHT'}", flush=True)
    return caught

def find(part, pred):
    for i, (nm, v) in enumerate(raw[part]):
        if pred(nm, v): return i
    raise KeyError

def nz(w): return sum(1 for a in w if a)
import re
def split(nm):
    m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]$', nm)
    return int(m.group(1)), m.group(2), [int(t) for t in m.group(3).split(',')]

res = []
res.append(run('1', raw['1'], "control: R5 unchanged"))
# 1 coefficient in R5 (a weight-6 term)
L = copy.deepcopy(raw['1']); i = find('1', lambda nm, v: nm.startswith('c^') and len(split(nm)[2]) == 6)
L[i][1] = str(Fr(L[i][1]) + Fr(1, 4)); res.append(run('1', L, f"M1 R5: coefficient of '{L[i][0]}' +1/4"))
# 2 coefficient in Psi3
L = copy.deepcopy(raw['z2']); i = find('z2', lambda nm, v: nm.startswith('c^') and len(split(nm)[2]) == 4)
L[i][1] = str(Fr(L[i][1]) * 2); res.append(run('z2', L, f"M2 Psi3: coefficient of '{L[i][0]}' doubled"))
# 2b constant 1/2 of Psi3 changed to 1/3 (only the low-order check can see it)
L = copy.deepcopy(raw['z2']); i = find('z2', lambda nm, v: nm == '1')
L[i][1] = '1/3'; res.append(run('z2', L, "M2b Psi3: constant 1/2 -> 1/3"))
# 3 coefficient in G
L = copy.deepcopy(raw['z2^2']); i = find('z2^2', lambda nm, v: nm == 'c^-3 He[0,1]')
L[i][1] = str(Fr(L[i][1]) + Fr(1, 32)); res.append(run('z2^2', L, f"M3 G: coefficient of 'c^-3 He[0,1]' -3/32 -> -2/32"))
# 4 dropped term in R5 (weight 3)
L = copy.deepcopy(raw['1']); i = find('1', lambda nm, v: nm.startswith('c^') and len(split(nm)[2]) == 3)
nm_d = L[i][0]; del L[i]; res.append(run('1', L, f"M4 R5: dropped term '{nm_d}'"))
# 5 exponent e changed (Psi3, a term with e = 0 -> e = -1)
L = copy.deepcopy(raw['z2']); i = find('z2', lambda nm, v: nm.startswith('c^0 He'))
e, k, w = split(L[i][0]); old = L[i][0]; L[i][0] = f"c^-1 He[{','.join(map(str, w))}]"
res.append(run('z2', L, f"M5 Psi3: '{old}' -> '{L[i][0]}'"))
# 5b exponent changed in R5 (e=1 -> e=-1, keeps parity)
L = copy.deepcopy(raw['1']); i = find('1', lambda nm, v: nm.startswith('c^1 He'))
e, k, w = split(L[i][0]); old = L[i][0]; L[i][0] = f"c^-1 He[{','.join(map(str, w))}]"
res.append(run('1', L, f"M5b R5: '{old}' -> '{L[i][0]}' (same parity)"))
# 6 swapped sign of the argument: He[w] -> He[-w] on a term with (-1)^(e+nz(w)) = -1 (else the pair is unchanged)
L = copy.deepcopy(raw['1']); i = find('1', lambda nm, v: nm.startswith('c^') and (split(nm)[0] + nz(split(nm)[2])) % 2 == 1)
e, k, w = split(L[i][0]); old = L[i][0]; L[i][0] = f"c^{e} He[{','.join(str(-a) for a in w)}]"
res.append(run('1', L, f"M6 R5: '{old}' -> '{L[i][0]}' (argument c -> -c inside H)"))
# 6b the same on G
L = copy.deepcopy(raw['z2^2']); i = find('z2^2', lambda nm, v: nm.startswith('c^') and ' ' in nm and (split(nm)[0] + nz(split(nm)[2])) % 2 == 1)
e, k, w = split(L[i][0]); old = L[i][0]; L[i][0] = f"c^{e} He[{','.join(str(-a) for a in w)}]"
res.append(run('z2^2', L, f"M6b G: '{old}' -> '{L[i][0]}'"))
# 7 an even pair replaced by its plain half c^e H_w(c) (breaks evenness; odd check)
L = copy.deepcopy(raw['z3']); i = find('z3', lambda nm, v: ' He[' in nm)
old = L[i][0]; L[i][0] = L[i][0].replace(' He[', ' H['); L.append([old.replace(' He[', ' H[').replace('[', '[-', 1).replace(',', ',-').replace('--', ''), "0"])
L = L[:-1]; res.append(run('z3', L, f"M7 [z3]: '{old}' -> plain half '{L[i][0]}'"))
# 8 ladder part: z5 coefficient
L = copy.deepcopy(raw['z5']); L[2][1] = str(Fr(L[2][1]) + Fr(1, 16)); res.append(run('z5', L, f"M8 [z5]: coefficient of '{L[2][0]}' +1/16"))
# 9 extra c^-2 constant in G (low orders only)
L = copy.deepcopy(raw['z2^2']); i = find('z2^2', lambda nm, v: nm == 'c^-2'); L[i][1] = '-1/2'
res.append(run('z2^2', L, "M9 G: constant c^-2 coefficient -3/4 -> -1/2"))
# 10 right part, wrong target: R5 against [z2] (sanity of the comparison)
terms = S.parse(raw['1']); d = lin((1, G.hpl_side(terms, 0)), (-1, P5[mon['z2']]))
print(f"{'M10 R5 compared with the [z2] target':70s} even-fails={bool(d)}")
print("control (unmutated R5) flagged by the checker:", res[0], "(must be False); mutants caught:", sum(res[1:]), "of", len(res) - 1)
