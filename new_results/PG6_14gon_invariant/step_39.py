import pickle
base = '<project>/research/08_theory/'
mons12 = pickle.load(open(base+'DOD1_dodecagon/mons.pkl','rb'))
print("mons12 size:", len(mons12))
# key format: (n_,(a,b),wc,letters) ; letters=(l2,l3,lu1) exponents
# total weight = a+b+sum(wc)+ (letter count). For 12gon letters are exponents of l2,l3,lu1
wd={}
for (n,ab,wc,let) in mons12:
    # letter degree: for 12gon letters stored as exponents; weight = a+b+sum(wc)+sum(let)?
    w = ab[0]+ab[1]+sum(wc)+sum(let)
    wd[w]=wd.get(w,0)+1
print("12gon weight dist (a+b+wc+sum(letters)):", wd)
# Show a few sample keys
for i,k in enumerate(list(mons12.keys())[:8]):
    print(k, "->", mons12[k])
# Check: are there degree<3 term-log monomials? look at E3 equivalent
# the max a+b+sum(wc) and the letters
maxw = max(ab[0]+ab[1]+sum(wc) for (n,ab,wc,let) in mons12)
print("max wlog degree:", maxw)
