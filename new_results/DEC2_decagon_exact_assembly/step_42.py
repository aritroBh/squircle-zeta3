import sympy as sp, pickle, os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
mine = pickle.load(open('E3n_sym_pq.pkl','rb'))
dec1 = pickle.load(open('<project>/research/08_theory/DEC1_decagon/E3n_exact_pq.pkl','rb'))
print("mine:", len(mine), " DEC1:", len(dec1))
mism = 0
only_mine = set(mine)-set(dec1); only_dec1 = set(dec1)-set(mine)
print("keys only in mine:", only_mine)
print("keys only in DEC1:", only_dec1)
worst = []
for k in dec1:
    if k not in mine: continue
    for j in range(4):
        a = mine[k][j]; b = dec1[k][j]
        da = sp.simplify(a[0]-b[0]); db = sp.simplify(a[1]-b[1])
        if da != 0 or db != 0:
            mism += 1
            worst.append((k,j,a,b))
print("coefficient mismatches:", mism)
for w in worst[:10]: print(w)
# DEC1 stores rows as list of 4 (p,q) with python rationals; mine sympy rationals -> compare as strings too
with open('progress.log','a') as f:
    f.write(f"step_03: exact E3n engine over Q(sqrt5): 55 monomials; monomial-by-monomial vs DEC1 E3n_exact_pq.pkl: {mism} mismatches\n")
print("PROGRESS: Task 1 comparison done, mismatches =", mism)
