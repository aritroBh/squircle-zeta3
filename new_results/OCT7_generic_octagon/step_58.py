import mpmath as mp, pickle, ast
mp.mp.dps = 700
d = pickle.load(open('lam_lll_450.pkl','rb'))
names = pickle.load(open('basis_names.pkl','rb'))['basis_names']
Lam = {k: mp.mpf(v) for k, v in d['Lam'].items()}

# reconstruct basis values at 700 digits
l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5)
z2 = mp.pi**2/6; z3 = mp.zeta(3)
logmap = {'l2': l2, 'l3': l3, 'l5': l5}
import re
def basis_value(name):
    # parse names like '1','l2','z2','Li2(-0.2)','Li3(0.1)l5','l2l5l3','z2l3'
    if name == '1': return mp.mpf(1)
    if name in logmap: return logmap[name]
    if name == 'z2': return z2
    if name == 'z3': return z3
    m = re.match(r'^(z2)?(Li[23]\((-?\d+\.\d+)\))?((l[235])*)$', name)
    val = mp.mpf(1)
    rest = name
    if rest.startswith('z2'): val *= z2; rest = rest[2:]
    mm = re.match(r'^Li([23])\((-?[\d.]+)\)', rest)
    if mm:
        order = int(mm.group(1)); arg = mp.mpf(mm.group(2))
        val *= mp.polylog(order, arg); rest = rest[mm.end():]
    for ch in re.findall(r'l[235]', rest):
        val *= logmap[ch]
    return val

for k, v in d['results'].items():
    if v[0] != 'REDUCES': continue
    rel = ast.literal_eval(v[1])
    nz = [(i,c) for i,c in enumerate(rel) if c != 0]
    tgt = Lam[k]
    res = mp.fsum([mp.mpf(c)*basis_value(names[i]) for i,c in nz if i < len(names)]) + sum(mp.mpf(c)*tgt for i,c in nz if i >= len(names))
    print(f"{k}: relation residual at 700 digits: {mp.nstr(abs(res), 3)}")
