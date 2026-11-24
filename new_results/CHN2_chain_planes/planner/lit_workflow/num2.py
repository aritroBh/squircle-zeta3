import mpmath as mp, sympy as sp, itertools, time
mp.mp.dps=60
U1,U2,U3,U4=sp.symbols('u1 u2 u3 u4')
def sym_key(T,B): return (tuple(sorted([sp.expand(x) for x in T],key=str)),tuple(sorted([sp.expand(x) for x in B],key=str)))
# BFS over Thomae with symbolic Gamma-prefactor tracking: prefactor = dict{expr:exp}
def thomae_step(T,B,i):
    d,e=B; a=T[i]; b,c=[T[j] for j in range(3) if j!=i]
    sg=sp.expand(d+e-a-b-c)
    fac={d:1,e:1,sg:1,a:-1,sp.expand(sg+b):-1,sp.expand(sg+c):-1}
    return ((sp.expand(d-a),sp.expand(e-a),sg),(sp.expand(sg+b),sp.expand(sg+c)),fac)
def mulfac(f,g):
    h=dict(f)
    for k,v in g.items(): h[k]=h.get(k,0)+v
    return {k:v for k,v in h.items() if v!=0}
start=((-U4,U3,1+U1+U3),(1+U3,1+U1+U2+U3))
target=sym_key((U2,U3,1+U1+U2+U3+U4),(1+U1+U2+U3,1+U2+U3+U4))
seen={sym_key(*start):(start[0],start[1],{})}; fr=[seen[sym_key(*start)]]; found=None
while fr and found is None:
    new=[]
    for T,B,f in fr:
        for i in range(3):
            T2,B2,g=thomae_step(T,B,i); k=sym_key(T2,B2)
            if k not in seen:
                seen[k]=(T2,B2,mulfac(f,g)); new.append(seen[k])
                if k==target: found=seen[k]
    fr=new
T2,B2,fac=found
print("Saalschutzian member:",T2,B2)
print("Thomae prefactor (3F2_orig = prod Gamma(arg)^exp * 3F2_member), before G:")
for k,v in fac.items(): print("   Gamma(%s)^%d"%(k,v))
def G(u1,u2,u3): return mp.gamma(1+u2)*mp.gamma(1+u1+u3)/mp.gamma(1+u1+u2+u3)
def hyp(T,B): return mp.hyp3f2(T[0],T[1],T[2],B[0],B[1],1)
u=[mp.mpf(x) for x in ('0.3','0.7','0.5','0.4')]
sub={U1:u[0],U2:u[1],U3:u[2],U4:u[3]}
t0=time.time()
orig=hyp([mp.mpf(sp.N(x.subs(sub),70)) for x in start[0]],[mp.mpf(sp.N(x.subs(sub),70)) for x in start[1]])
mem=hyp([mp.mpf(sp.N(x.subs(sub),70)) for x in T2],[mp.mpf(sp.N(x.subs(sub),70)) for x in B2])
pref=mp.mpf(1)
for k,v in fac.items(): pref*=mp.gamma(mp.mpf(sp.N(k.subs(sub),70)))**v
print("generic point: |orig - pref*member| =",mp.nstr(abs(orig-pref*mem),5)," (orig=%s)"%mp.nstr(orig,20)," time %.1fs"%(time.time()-t0))
# Vol/16 on plane Pi1: u=(al,al,be,be) -> compare with Gamma-factor * F(al,be)
al,be=mp.mpf('0.35'),mp.mpf('0.6')
V=G(al,al,be)*hyp([-be,be,1+al+be],[1+be,1+2*al+be])
Fst=hyp([al,be,1+2*al+2*be],[1+2*al+be,1+al+2*be])
subp={U1:al,U2:al,U3:be,U4:be}
prefp=G(al,al,be)
for k,v in fac.items(): prefp*=mp.gamma(mp.mpf(sp.N(k.subs(subp),70)))**v
print("Pi1: |Vol/16 - Gamma-factor*F(al,be)| =",mp.nstr(abs(V-prefp*Fst),5),"  Vol/16=",mp.nstr(V,25))
al,be=mp.mpf('0.35'),mp.mpf('0.6')
V2=G(al,be,al)*hyp([-be,al,1+2*al],[1+al,1+al+be+al])
subp2={U1:al,U2:be,U3:al,U4:be}
prefp2=G(al,be,al)
for k,v in fac.items(): prefp2*=mp.gamma(mp.mpf(sp.N(k.subs(subp2),70)))**v
print("Pi2: |Vol/16 - Gamma-factor*F(be,al)| =",mp.nstr(abs(V2-prefp2*hyp([be,al,1+2*al+2*be],[1+2*be+al,1+be+2*al])),5))
# Nonterminating Saalschutz (Gasper-Rahman II.24 classical limit, as in Koornwinder (4.1)) for F(s,t)
def nts(a,b,c,e):
    f=a+b+c-e+1
    lhs1=mp.hyp3f2(a,b,c,e,f,1)
    pref2=mp.gamma(e-1)*mp.gamma(a-e+1)*mp.gamma(b-e+1)*mp.gamma(c-e+1)*mp.gamma(f)/(mp.gamma(1-e)*mp.gamma(a)*mp.gamma(b)*mp.gamma(c)*mp.gamma(a+b+c-2*e+2))
    part=mp.hyp3f2(a-e+1,b-e+1,c-e+1,2-e,a+b+c-2*e+2,1)
    rhs=mp.gamma(a-e+1)*mp.gamma(b-e+1)*mp.gamma(c-e+1)*mp.gamma(f)/(mp.gamma(1-e)*mp.gamma(b+c-e+1)*mp.gamma(a+c-e+1)*mp.gamma(a+b-e+1))
    return lhs1,pref2,part,rhs
for (s,t) in [(mp.mpf('0.35'),mp.mpf('0.6')),(mp.mpf('0.2'),mp.mpf('0.2')),(mp.mpf('0.15'),mp.mpf('0.45'))]:
    for which,e in ((0,1+2*s+t),(1,1+s+2*t)):
        lhs1,pref2,part,rhs=nts(s,t,1+2*s+2*t,e)
        print("(s,t)=(%s,%s) e=bottom[%d]: |F + B*F2 - A| = %s ; F=%s ; B*F2=%s ; A=%s ; F2(partner)=%s"%(s,t,which,mp.nstr(abs(lhs1+pref2*part-rhs),5),mp.nstr(lhs1,15),mp.nstr(pref2*part,15),mp.nstr(rhs,15),mp.nstr(part,15)))
# ratio of partner to F(s,t) at two points to see if it is a Gamma ratio: just report
