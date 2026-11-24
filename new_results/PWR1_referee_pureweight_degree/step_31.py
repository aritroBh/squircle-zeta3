import pwr1_engine as E, sympy as sp, re

txt = open('../../10_frontier/H1_octahedron_end/logform_general_n.out').read()
# parse "[eta^k] log B_n = <expr>   (= ...)" lines per n-block
blocks = re.split(r'n=(\d+) \(K=12\)', txt)
# blocks[0] empty, then pairs (n, text)
zvars = {f'z{k}': E.Z[k] for k in range(2, 13)}
ok_all = True
for i in range(1, len(blocks), 2):
    n = int(blocks[i]); body = blocks[i+1]
    if n not in (2, 3): continue
    lB, B = E.logB(n)
    for k in range(1, 9):
        m = re.search(rf'\[eta\^{k}\] log B_n = (.*?)(?:   =|\n)', body)
        theirs = sp.sympify(m.group(1).strip(), locals=zvars)
        mine = lB.coeff(k)
        same = sp.expand(mine - theirs) == 0
        if not same: ok_all = False
        print(f"n={n} eta^{k}: exact match {same}")
print("ALL logB match:", ok_all)
open('progress.log','a').write("step_18: n=3..6 engine: gamma-free, weight-homogeneous, w1=0; w2/w3 match H1 T6 A2/A3; logB(n=2,3) EXACT match vs H1 logform_general_n.out through eta^8\n")
