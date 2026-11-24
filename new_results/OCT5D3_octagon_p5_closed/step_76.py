import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
src = open("oct5d3_basis.py").read()
addition = '''

def atom_specs():
    """(label, kind, spec) with spec = list of (base_atom_name, power) for numerical evaluation.
    Base names match oct5d3_atoms.base_atoms(): lam, L, pi2, z2..z5, Lchi3, Lchi5,
    Li2(b), Li3(b/-b/1/sqrt2), Li4(b/-b/-b^2/1/sqrt2/1/2/1-1/sqrt2/2-sqrt2),
    Li5(14 args), Li32(b,1), Li32(1/sqrt2,1)."""
    specs = []
    for name, z in LI5_ARGS:
        specs.append((f"Li5({name})", 'vis', [(f"Li5({name})", 1)]))
    for i in range(5, -1, -1):
        lab = f"lam^{i} L^{5-i}" if i and 5-i else ("lam^5" if i == 5 else "L^5")
        specs.append((lab, 'vis', ([('lam', i)] if i else []) + ([('L', 5-i)] if 5-i else [])))
    for i in range(3, -1, -1):
        specs.append((f"lam^{i} L^{3-i} pi2", 'inv',
                      ([('lam', i)] if i else []) + ([('L', 3-i)] if 3-i else []) + [('pi2', 1)]))
    for i in range(3, -1, -1):
        specs.append((f"lam^{i} L^{3-i} Li2(b)", 'vis',
                      ([('lam', i)] if i else []) + ([('L', 3-i)] if 3-i else []) + [('Li2(b)', 1)]))
    specs.append(("lam pi2^2", 'inv', [('lam', 1), ('pi2', 2)]))
    specs.append(("L pi2^2", 'inv', [('L', 1), ('pi2', 2)]))
    specs.append(("lam pi2 Li2(b)", 'inv', [('lam', 1), ('pi2', 1), ('Li2(b)', 1)]))
    specs.append(("L pi2 Li2(b)", 'inv', [('L', 1), ('pi2', 1), ('Li2(b)', 1)]))
    specs.append(("lam Li2(b)^2", 'vis', [('lam', 1), ('Li2(b)', 2)]))
    specs.append(("L Li2(b)^2", 'vis', [('L', 1), ('Li2(b)', 2)]))
    for i in range(2, -1, -1):
        specs.append((f"lam^{i} L^{2-i} z3", 'inv',
                      ([('lam', i)] if i else []) + ([('L', 2-i)] if 2-i else []) + [('z3', 1)]))
    for i in range(2, -1, -1):
        for name, z in LI3_ARGS:
            specs.append((f"lam^{i} L^{2-i} Li3({name})", 'vis',
                          ([('lam', i)] if i else []) + ([('L', 2-i)] if 2-i else []) + [(f"Li3({name})", 1)]))
    specs.append(("pi2 z3", 'inv', [('pi2', 1), ('z3', 1)]))
    for name, z in LI3_ARGS:
        specs.append((f"pi2 Li3({name})", 'inv', [('pi2', 1), (f"Li3({name})", 1)]))
    specs.append(("Li2(b) z3", 'inv', [('Li2(b)', 1), ('z3', 1)]))
    for name, z in LI3_ARGS:
        specs.append((f"Li2(b) Li3({name})", 'vis', [('Li2(b)', 1), (f"Li3({name})", 1)]))
    for name, z in LI4_ARGS:
        specs.append((f"lam Li4({name})", 'vis', [('lam', 1), (f"Li4({name})", 1)]))
        specs.append((f"L Li4({name})", 'vis', [('L', 1), (f"Li4({name})", 1)]))
    specs.append(("z5", 'inv', [('z5', 1)]))
    specs.append(("Lchi8_5", 'inv', [('Lchi5', 1)]))
    return specs

C_SPECS = [("Li32(b,1)", [("Li32(b,1)", 1)], C_TERMS[0][2]),
           ("Li32(1/sqrt2,1)", [("Li32(1/sqrt2,1)", 1)], C_TERMS[1][2])]
'''
open("oct5d3_basis.py","a").write(addition)
import subprocess, sys
r = subprocess.run([sys.executable, "-c", """
import sys; sys.path.insert(0,'.')
from oct5d3_basis import atom_specs, atom_list
sp = atom_specs(); al = atom_list()
print(len(sp), len(al))
assert [s[0] for s in sp] == [a[0] for a in al], 'label mismatch'
assert [s[1] for s in sp] == [a[1] for a in al], 'kind mismatch'
print('labels+kinds match')
"""], capture_output=True, text=True)
print(r.stdout, r.stderr[-800:] if r.returncode else "")
