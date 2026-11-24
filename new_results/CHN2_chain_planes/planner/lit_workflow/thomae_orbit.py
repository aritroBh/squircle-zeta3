import sympy as sp, itertools
u1,u2,u3,u4,al,be=sp.symbols('u1 u2 u3 u4 al be')
top=(-u4,u3,1+u1+u3); bot=(1+u3,1+u1+u2+u3)
def key(T,B):
    return (tuple(sorted([sp.expand(x) for x in T],key=str)),tuple(sorted([sp.expand(x) for x in B],key=str)))
def thomae(T,B):
    d,e=B; out=[]
    for i in range(3):
        a=T[i]; b,c=[T[j] for j in range(3) if j!=i]
        sig=sp.expand(d+e-a-b-c)
        out.append(((d-a,e-a,sig),(sig+b,sig+c)))
    return out
seen={key(top,bot):(top,bot)}
frontier=[(top,bot)]
while frontier:
    new=[]
    for T,B in frontier:
        for T2,B2 in thomae(T,B):
            k=key(T2,B2)
            if k not in seen:
                seen[k]=(T2,B2); new.append((T2,B2))
    frontier=new
print("orbit size",len(seen))
def is_saal(T,B,sub=None):
    res=[]
    for i,j in itertools.permutations(range(3),2):
        k=3-i-j
        p,q,r=T[i],T[j],T[k]
        for B0 in (B,B[::-1]):
            ok1=sp.simplify(r-(1+2*p+2*q))==0
            ok2=sp.simplify(B0[0]-(1+2*p+q))==0 and sp.simplify(B0[1]-(1+p+2*q))==0
            if ok1 and ok2: res.append((p,q))
    return res
planes={'Pi1':{u1:al,u2:al,u3:be,u4:be},'Pi2':{u1:al,u3:al,u2:be,u4:be}}
for name,sub in planes.items():
    print("==",name)
    for k,(T,B) in seen.items():
        T2=[sp.expand(x.subs(sub)) for x in T]; B2=[sp.expand(x.subs(sub)) for x in B]
        r=is_saal(T2,B2)
        print(T2,B2,"SAAL(p,q)=",r)
print("== general (no plane)")
for k,(T,B) in seen.items():
    print(T,B,is_saal(T,B))
