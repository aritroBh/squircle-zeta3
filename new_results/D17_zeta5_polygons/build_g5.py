
import sympy as sp, sys, time, pickle
sys.path.insert(0, "research/08_theory/D17_zeta5_polygons")
import sym_D5
t0=time.time()
D, syms = sym_D5.build_D(5, reduced=True)
print("build time %.1fs"%(time.time()-t0))
L0,L1 = syms['L0'], syms['L1']
D5 = sp.expand(D[5])
print("D5 terms:", len(D5.as_ordered_terms()))
t0=time.time()
G = {}
for (a,b) in [(4,1),(3,2),(2,3),(1,4),(5,0),(0,5)]:
    g = sp.expand(D5.coeff(L0,a).coeff(L1,b))
    G[(a,b)] = g
    print("G%d%d: %d terms  (%.1fs)"%(a,b,len(g.as_ordered_terms()),time.time()-t0))
out = {f"{a}{b}": sp.srepr(G[(a,b)]) for (a,b) in G}
with open("G5_reduced.txt","w") as f:
    for k,v in out.items():
        f.write(k+"\n"+v+"\n")
print("saved G5_reduced.txt")
with open("D0B.txt","w") as f:
    f.write(sp.srepr(syms['D0'])+"\n###\n"+sp.srepr(syms['B']))
print("saved D0B.txt")
