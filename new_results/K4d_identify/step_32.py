import sympy as sp, itertools
rt2 = sp.sqrt(2)
val4=[]
for ax in range(3):
    for s in [1,-1]:
        v=[0,0,0]; v[ax]=s*rt2/4; val4.append(tuple(v))
val6=[tuple(s*rt2/6 for s in signs) for signs in itertools.product([1,-1],repeat=3)]
verts=val4+val6
NA=[(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]
NA=[tuple(sp.Rational(x)/rt2 for x in n) for n in NA]
def avec(v): return [sp.simplify(sum(n[i]*v[i] for i in range(3))) for n in NA]
# representative apex-base edge: val4[1]=(rt2/4,0,0) and find a val6 neighbor at dist^2 1/8
def dist2(a,b): return sp.simplify(sum((a[i]-b[i])**2 for i in range(3)))
print("=== APEX-BASE edges (first vertex valence-4, term values) ===")
for va in val4:
    for vb in val6:
        if dist2(va,vb)==sp.Rational(1,8):
            aa=avec(va); ab=avec(vb)
            print("Va=",va," a=",aa)
            print("Vb=",vb," a=",ab)
            # ratio for terms nonzero at both
            for k in range(6):
                if aa[k]!=0 and ab[k]!=0:
                    print("   term %d nonzero both: %s -> %s  ratio=%s"%(k,aa[k],ab[k],sp.simplify(ab[k]/aa[k])))
            print()
    break  # just first valence-4
