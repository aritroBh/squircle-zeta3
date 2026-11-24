
import pickle, re
fs = pickle.load(open('../COR8_depth_two_any_angle/full_sols_v3.pkl','rb'))
alln = fs['alln']
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
bad = set()
for item in alln:
    if item[0]=='Li4': continue
    s = item[1]
    parts = []
    for fac in split_factors(s):
        fac = fac.strip()
        if fac.endswith('**2'): parts += [fac[:-3]]*2
        else: parts.append(fac)
    for fac in parts:
        if not re.match(r'(Li2|Li3|log)\((.*)\)$', fac):
            bad.add(fac)
print("unparseable factors:", bad)
# also distinct function names
names = set()
for item in alln:
    if item[0]=='Li4': continue
    for fac in split_factors(item[1]):
        fac = fac.strip().rstrip('*2') if fac.strip().endswith('**2') else fac.strip()
        m = re.match(r'(\w+)\(', fac)
        if m: names.add(m.group(1))
print("function names in prods:", names)
