import pickle
wd = '<project>/research/08_theory/K1f_hexagon_order5'
polys = pickle.load(open(f'{wd}/E5_polys.pkl','rb'))
E5 = polys[5]
abc = sorted({(a,b,c) for (n,a,b,c) in E5})
nl_abc = sorted({(a,b,c) for (n,a,b,c) in E5 if c>=1})
loc_abc = sorted({(a,b) for (n,a,b,c) in E5 if c==0})
print("E[5] monomials:", len(E5))
print("max t-degree:", max(n for (n,a,b,c) in E5))
print("distinct (a,b,c):", len(abc))
print("(a,b,c) with c>=1:", nl_abc)
print("(a,b) local:", loc_abc)
nmax_per_abc = {}
for (n,a,b,c) in E5:
    nmax_per_abc[(a,b,c)] = max(nmax_per_abc.get((a,b,c),-1), n)
print("n-range per (a,b,c) c>=1:")
for k in nl_abc:
    print("  ", k, "n=0..%d"%nmax_per_abc[k])
# also lower orders for cross-checks
for j in range(5):
    pj = polys[j]
    print(f"E[{j}]: {len(pj)} monomials, max t:", max(n for (n,a,b,c) in pj), "abc:", sorted({(a,b,c) for (n,a,b,c) in pj})[:12])
