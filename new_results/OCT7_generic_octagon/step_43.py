import sympy as sp, pickle
Jtab = pickle.load(open("research/08_theory/K1d_engine_fixed/Jtab_t.pkl","rb"))
logs_found = set()
for k,v in Jtab.items():
    for atom in sp.preorder_traversal(v):
        if isinstance(atom, sp.log):
            logs_found.add(atom)
print("log atoms in Jtab:", logs_found)
Jloc = pickle.load(open("research/08_theory/OCT5_octagon_exact/Jloc.pkl","rb"))
logs_found2 = set()
for k,v in Jloc.items():
    for atom in sp.preorder_traversal(v):
        if isinstance(atom, sp.log):
            logs_found2.add(atom)
print("log atoms in Jloc:", logs_found2)
