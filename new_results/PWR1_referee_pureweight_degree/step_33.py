import pwr1_engine as E, sympy as sp, re, json
txt = open('../../10_frontier/H1_octahedron_end/logform_general_n.out').read()
blocks = re.split(r'n=(\d+) \(K=12\)', txt)
zvars = {f'z{k}': E.Z.get(k, sp.Symbol(f'z{k}')) for k in range(2, 13)}
zvars['pi'] = sp.pi
ok_all = True
for i in range(1, len(blocks), 2):
    n = int(blocks[i]); body = blocks[i+1]
    if n not in (2, 3): continue
    lB, B = E.logB(n)
    for k in range(2, 9):
        m = re.search(rf'\[eta\^{k}\] log B_n = (.*?)(?:   =|\n)', body)
        if m is None: print(f"n={n} k={k}: line not found"); ok_all=False; continue
        theirs = sp.sympify(m.group(1).strip(), locals=zvars)
        same = sp.expand(lB.coeff(k) - theirs) == 0
        ok_all &= same
        if not same: print(f"  MISMATCH n={n} eta^{k}")
print("logB(n=2,3) exact match vs H1 through eta^8:", ok_all)

j3 = json.load(open('../../10_frontier/H1_octahedron_end/routeA_n3_K12.json'))
print(list(j3.keys()), list(j3['3'].keys()))
print(str(j3['3']['eta'][:4]))
