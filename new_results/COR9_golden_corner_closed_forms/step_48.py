import pickle, sympy as sp
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
d = pickle.load(open(WD+"golden_B4.pkl", "rb"))
B4_NAMES = d["names"]; nb = len(B4_NAMES)
pure = pickle.load(open(WD+"golden_B4_relations.pkl", "rb"))
# rational relations: rows with sqrt5-part zero
rats = []
for g in pure:
    if all(c == 0 for c in g[nb:]):
        rats.append(g[:nb])
    elif all(c == 0 for c in g[:nb]):
        rats.append(g[nb:])
R = sp.Matrix(rats)
print("rational relations:", R.rows, " rank:", R.rank())
pivots = R.rref()[1]
print("pivot (dependent) columns:", [B4_NAMES[i] for i in pivots])
free = [i for i in range(nb) if i not in pivots]
print(f"independent basis ({len(free)}):")
for i in free: print("   ", B4_NAMES[i])
