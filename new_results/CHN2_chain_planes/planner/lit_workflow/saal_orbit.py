import sympy as sp, itertools
s,t=sp.symbols('s t')
def key(T,B): return (tuple(sorted([sp.expand(x) for x in T],key=str)),tuple(sorted([sp.expand(x) for x in B],key=str)))
def thomae(T,B):
    d,e=B; out=[]
    for i in range(3):
        a=T[i]; b,c=[T[j] for j in range(3) if j!=i]
        sg=sp.expand(d+e-a-b-c); out.append(((d-a,e-a,sg),(sg+b,sg+c)))
    return out
def orbit(T,B):
    seen={key(T,B):(T,B)}; fr=[(T,B)]
    while fr:
        new=[]
        for T1,B1 in fr:
            for T2,B2 in thomae(T1,B1):
                k=key(T2,B2)
                if k not in seen: seen[k]=(T2,B2); new.append((T2,B2))
        fr=new
    return list(seen.values())
def saal_family(T,B):
    res=[]
    for i,j in itertools.permutations(range(3),2):
        k=3-i-j; p,q,r=T[i],T[j],T[k]
        for B0 in (B,B[::-1]):
            if sp.expand(r-(1+2*p+2*q))==0 and sp.expand(B0[0]-(1+2*p+q))==0 and sp.expand(B0[1]-(1+p+2*q))==0: res.append((p,q))
    return res
F=((s,t,1+2*s+2*t),(1+2*s+t,1+s+2*t))
# nonterminating Saalschutz partner for e = 1+2s+t  and e = 1+s+2t
def partner(T,B,which):
    a,b,c=T; e=B[which]; f=B[1-which]
    return ((sp.expand(a-e+1),sp.expand(b-e+1),sp.expand(c-e+1)),(sp.expand(2-e),sp.expand(f-e+1)))
for which in (0,1):
    P=partner(*F,which)
    print("partner (e=bottom[%d]):"%which,P)
    ob=orbit(*P)
    print("  orbit size",len(ob))
    for T2,B2 in ob:
        print("   ",T2,B2,"saal-family match:",saal_family(T2,B2))
print("Thomae orbit of F itself:")
for T2,B2 in orbit(*F): print("   ",T2,B2,saal_family(T2,B2))
