import sympy as sp, itertools
rt2 = sp.sqrt(2)
val4=[]
for ax in range(3):
    for s in [1,-1]:
        v=[0,0,0]; v[ax]=s*rt2/4; val4.append(tuple(v))
val6=[tuple(s*rt2/6 for s in signs) for signs in itertools.product([1,-1],repeat=3)]
NA=[(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]
NA=[tuple(sp.Rational(x)/rt2 for x in n) for n in NA]
def avec(v): return [sp.simplify(sum(n[i]*v[i] for i in range(3))) for n in NA]
def dist2(a,b): return sp.simplify(sum((a[i]-b[i])**2 for i in range(3)))
print("=== BASE-BASE edges (valence-6 to valence-6, dist^2=2/9) ===")
cnt=0
for i in range(8):
    for j in range(i+1,8):
        va,vb=val6[i],val6[j]
        if dist2(va,vb)==sp.Rational(2,9):
            aa,ab=avec(va),avec(vb)
            print("Va=",va," a=",aa)
            print("Vb=",vb," a=",ab)
            for k in range(6):
                if aa[k]!=0 and ab[k]!=0:
                    print("   term %d both: %s -> %s  ratio=%s"%(k,aa[k],ab[k],sp.simplify(ab[k]/aa[k])))
            print()
            cnt+=1
    if cnt>=3: break
