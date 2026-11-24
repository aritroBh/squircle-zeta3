WD='<project>/research/08_theory/FLD1_field_rule_proof'
with open(WD+'/progress.log','a') as fh:
    fh.write("step28: weight bound wt(z)+|w|=k+1 SATURATED, 0 violations; letters={0,1,-1} (-1 needed ~87%); denominators only powers of c\n")
    fh.write("step29: k=4 construction == G2_M4alpha B4 (R4,Psi2,z3=-B3z2,z4=F1,z2^2=F2) as functions. Consistent.\n")
import pickle
Z = pickle.load(open(WD+'/mu_zsum_K5.pkl','rb'))
def wt(zkey): return sum((i+2)*e for i,e in enumerate(zkey))
print("Z-sum weight bound check: max |w| <= k - wt(z) ?")
for k in range(6):
    for zkey, terms in sorted(Z[k].items()):
        maxw = max((sum(abs(x) for x in w) for w,c in terms), default=0)
        # also pole orders
        import sympy as sp
        n=sp.Symbol('n')
        print(f'  k={k} z={zkey} (wt {wt(zkey)}): max Z-sum weight = {maxw},  k-wt(z) = {k-wt(zkey)},  ok={maxw <= k-wt(zkey)}')