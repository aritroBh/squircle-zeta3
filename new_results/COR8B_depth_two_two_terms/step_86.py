
import pickle, re, math, json
import mpmath as mp
import sympy as sp
mp.mp.dps = 40
c_sym = sp.symbols('c')
fs = pickle.load(open('../COR8_depth_two_any_angle/full_sols_v3.pkl','rb'))
alln, sols = fs['alln'], fs['sols']
def split_factors(s):
    parts, depth, cur = [], 0, ''
    for ch in s:
        if ch == '(': depth += 1; cur += ch
        elif ch == ')': depth -= 1; cur += ch
        elif ch == '*' and depth == 0:
            if cur: parts.append(cur); cur = ''
        else: cur += ch
    if cur: parts.append(cur)
    return parts
def ev_factor(fac, c0, c0s):
    m = re.match(r'(Li2|Li3|Li4|log)\((.*)\)$', fac)
    fname, arg = m.group(1), m.group(2)
    z = sp.N(sp.sympify(arg).subs(c_sym, sp.Rational(c0s)), 40)
    z = mp.mpc(z)
    if fname=='Li2': return mp.polylog(2, z)
    if fname=='Li3': return mp.polylog(3, z)
    if fname=='Li4': return mp.polylog(4, z)
    if fname=='log': return mp.log(z)
def ev_item(item, c0, c0s):
    if item[0]=='Li4':
        return ev_factor('Li4(%s)'%item[1], c0, c0s)
    s = item[1].strip()
    dup = False
    if s.endswith('**2'): s = s[:-3].strip(); dup = True
    parts = [f.strip() for f in split_factors(s)]
    if dup: parts.append(parts[-1])
    p = mp.mpc(1)
    for fac in parts:
        p *= ev_factor(fac, c0, c0s)
    return p
def mpf_f(fr): return mp.mpf(fr.numerator)/mp.mpf(fr.denominator)
for c0s in ['1/2']:
    c0 = mp.mpf(c0s)
    tot = mp.mpc(0)
    for e in (-1,0,1):
        for i, item in enumerate(alln):
            cf = sols[e][i]
            if not cf: continue
            tot += (c0**e) * mpf_f(cf) * ev_item(item, c0, c0s)
    print("COR8 classical at c=1/2:", mp.nstr(tot, 30))
print("Omega - DT_A at c=1/2 =", -0.0416576932572766672882308026250919127289839427 + 0.370468167158478130188006668187879050052512399)
