verify = r'''"""verify_pg6.py - PG6: 14-gon P3(1) resolution and the (2,1) invariant delta(N3).
Recomputes from saved exact data (geom14/E3n14/letters14/mons14) + engine14/sym14. Run: PY verify_pg6.py
"""
import sys, pickle, itertools, mpmath as mp
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import engine14 as E
import sym14 as S
import sympy as sp
mp.mp.dps=50
PASS=[True]
def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")
    PASS[0]=PASS[0] and bool(cond)
g = S.init()
m=7
def fsum(th,p): return sum(mp.exp(p*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
def r_rp(th,p):
    f=fsum(th,p); fp=sum(-mp.sin(th-i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    return f**(-1/p), -(1/p)*f**(-1/p-1)*fp
def ig(th,p):
    r,rp=r_rp(th,p); return mp.sqrt(r*r+rp*rp)
def perim(p, chamber):
    return 2*m*mp.quad(lambda th: ig(th,p), chamber, method='tanh-sinh', maxdegree=8)
c3={-4:-mp.mpf(7)/240,-3:mp.mpf(3)/10,-2:-mp.mpf(169)/120,-1:mp.mpf(61)/30,1:-mp.mpf(61)/30,2:mp.mpf(169)/120,3:-mp.mpf(3)/10,4:mp.mpf(7)/240}
def P3_fd(chamber,h):
    vals={k:perim(1+k*h,chamber) for k in range(-4,5)}
    return sum(c3[k]*vals[k] for k in c3)/h**3
TRUTH=mp.mpf('-5.5080836980534783'); GOLD=mp.mpf('-5.5080732424212877')
p_fd = P3_fd([-mp.pi/14,mp.pi/14], mp.mpf('0.0025'))
p_gold = P3_fd([mp.pi/2+mp.mpf('1e-9'), mp.pi/2+mp.pi/m-mp.mpf('1e-9')], mp.mpf('0.0025'))
E3=pickle.load(open('E3n14.pkl','rb'))['E3n']; ORDER=[3,4,1,2,5,6,0]
apn=[[E.fev(g['data'][i][0],50),E.fev(g['data'][i][1],50)] for i in ORDER]
LL=mp.sqrt(E.fev(g['L'],50))
def E3ev(t):
    logs=[mp.log(apn[i][0]+apn[i][1]*t) for i in range(7)]; v=mp.mpf(0)
    for k,c in E3.items():
        poly=((E.fev(c[3],50)*t+E.fev(c[2],50))*t+E.fev(c[1],50))*t+E.fev(c[0],50)
        lv=mp.mpf(1)
        for i,e in enumerate(k):
            if e: lv*=logs[i]**e
        v+=poly*lv
    return v
Ie=mp.quad(E3ev,[0,mp.mpf('0.5'),1],method='tanh-sinh',maxdegree=8)
p_chain=84*LL*Ie
check("A1: centered-chamber FD ~ truth", abs(p_fd-TRUTH)<1e-8, f"fd={mp.nstr(p_fd,16)}")
check("A2: 84*L*int E3n = truth", abs(p_chain-TRUTH)<1e-12, f"chain={mp.nstr(p_chain,20)}")
check("A3: gold chamber reproduces WRONG gold (bug demo)", abs(p_gold-GOLD)<1e-8 and abs(p_gold-TRUTH)>1e-6, f"goldchamber={mp.nstr(p_gold,16)}")
mons14=pickle.load(open('mons14.pkl','rb'))
INTC={}
def INTc(n,a,b,wc):
    k=(n,a,b,wc)
    if k not in INTC: INTC[k]=E.INT(n,a,b,wc)
    return INTC[k]
S14={}; z3loc=E.ZERO
for (n,ab,wc,fk),coef in mons14.items():
    J=INTc(n,ab[0],ab[1],wc); nfam=sum(1 for c in wc if c>0)
    lm={}
    if fk[0]: lm['l2']=fk[0]
    if fk[1]: lm['la']=fk[1]
    if fk[2]: lm['lu']=fk[2]
    lm=E.mkey(lm)
    for mon,c in J.items():
        S14[E.mdot(mon,lm)]=E.fadd(S14.get(E.mdot(mon,lm),E.ZERO), E.fmul(c,coef))
        d=dict(mon)
        if nfam==0 and d.get('z3')==1 and all(s in('z3','l2','la','lu') for s in d):
            z3loc=E.fadd(z3loc, E.fmul(c,coef))
l2v=mp.log(2); lav=mp.log(2*mp.cos(2*mp.pi/7)); luv=mp.log(2*mp.cos(mp.pi/7)); z2v=mp.pi**2/6; z3v=mp.zeta(3)
Xv={i:E.fev(E.X[i],50) for i in E.IDX}
w={i:(lambda tt,i=i: mp.log(1-Xv[i]*tt)) for i in E.IDX}; u={i:(lambda tt,i=i: 1-Xv[i]*tt) for i in E.IDX}
def Li21(x,y):
    tot=mp.mpf(0);ps=mp.mpf(0);xp=mp.mpf(1);yp=mp.mpf(1)
    for mm in range(1,4000):
        xp*=x
        if mm>=2: tot+=xp*ps/mm**2
        yp*=y; ps+=yp/mm
    return tot
def Aser(x):
    tot=mp.mpf(0);xp=mp.mpf(1);H=mp.mpf(0)
    for n in range(1,6000): xp*=x;H+=mp.mpf(1)/n;tot+=H*xp/n**2
    return tot
NUM={'l2':l2v,'la':lav,'lu':luv,'z2':z2v,'z3':z3v}
for i in E.IDX: NUM['L2x%d'%i]=mp.polylog(2,Xv[i]);NUM['L3x%d'%i]=mp.polylog(3,Xv[i]);NUM['Ax%d'%i]=Aser(Xv[i])
for (i,j) in E.Lam: NUM['Lam%d_%d'%(i,j)]=Li21(Xv[i],Xv[j]/Xv[i])
for (i,j) in E.T2s: NUM['T2_%d_%d'%(i,j)]=mp.quad(lambda tt:w[j](tt)/u[i](tt),[0,1])
for (i,j) in E.T3ps: NUM['T3p_%d_%d'%(i,j)]=mp.quad(lambda tt:w[j](tt)**2/u[i](tt),[0,1])
for (i,j,k) in E.T3a: NUM['T3_%d_%d%d'%(i,j,k)]=mp.quad(lambda tt:w[j](tt)*w[k](tt)/u[i](tt),[0,1])
def e_eval(P):
    tot=mp.mpf(0)
    for mon,c in P.items():
        v=E.fev(c,50)
        for (s,ex) in mon: v*=NUM[s]**ex
        tot+=v
    return tot
Snum=e_eval(S14); p_asm=84*LL*Snum
check("B1: 84*L*S(engine assembly) = truth", abs(p_asm-TRUTH)<1e-35, f"diff={mp.nstr(abs(p_asm-TRUTH),3)}")
C7=48*7*mp.sin(mp.pi/14)**4*mp.cos(mp.pi/14)*(3*mp.sin(mp.pi/14)**2-1)
check("B2: 84*L*zeta3_local = C_7", abs(84*LL*E.fev(z3loc,50)-C7)<1e-40, f"diff={mp.nstr(abs(84*LL*E.fev(z3loc,50)-C7),3)}")
def is_w3(nm): return nm.startswith('L3x') or nm.startswith('Ax') or nm.startswith('Lam') or nm.startswith('T3p') or nm.startswith('T3_')
def atom_sym(name):
    if name.startswith('L3x'): return S.sym_L3x(int(name[3:]))
    if name.startswith('Ax'): return S.sym_Ax(int(name[2:]))
    if name.startswith('Lam'): return S.sym_Lam(int(name[3]),int(name[5]))
    if name.startswith('T3p'): return S.sym_T3p(int(name[4]),int(name[6]))
    if name.startswith('T3_'): return S.sym_T3a(int(name[3]),int(name[5]),int(name[6]))
sigma={}
for mon,c in S14.items():
    syms=[s for s in dict(mon) if s not in ('l2','la','lu','z2','z3')]
    others=[s for s in dict(mon) if s in ('l2','la','lu','z2','z3')]
    if len(syms)==1 and not others and is_w3(syms[0]):
        sym,scalar=atom_sym(syms[0]); tot=E.fmul(c,scalar)
        for key,cf in sym.items(): sigma[key]=E.fadd(sigma.get(key,E.ZERO), E.fs(tot,cf))
letters=[0,1,2]; basis3=list(itertools.product(letters,repeat=3)); bidx={t:i for i,t in enumerate(basis3)}
Dg=[]
for a in letters:
    for b in letters:
        for c in letters:
            v={}
            for t in [(a,b,c),(b,a,c),(b,c,a)]: v[t]=v.get(t,0)+1
            Dg.append(v)
def t2v(T):
    v=[Fr(0)]*27
    for t,c in T.items():
        if t in bidx: v[bidx[t]]+=Fr(c)
    return v
L=pickle.load(open('letters14.pkl','rb')); CRs=[e[0] for e in L['CRset']]
LI3v=[t2v(S.li3_symbol(z)) for z in CRs]; Dgv=[t2v(T) for T in Dg]
M=sp.Matrix(27,len(LI3v)+len(Dgv),lambda i,j: sp.Rational((LI3v+Dgv)[j][i]))
rankM=M.rank()
red_ok=True; sols=[]
for k in range(3):
    sc=[Fr(0)]*27
    for key,cFT in sigma.items():
        A,B,C,d=cFT
        if key in bidx: sc[bidx[key]]+= [Fr(A)/d,Fr(B)/d,Fr(C)/d][k]
    tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(sc[i]))
    if M.row_join(tgt).rank()!=rankM: red_ok=False
    sols.append(M.gauss_jordan_solve(tgt)[0])
check("B3: sigma(S) reduces to 24-CR Li3 basis mod Decomp (alphabet rule)", red_ok, f"rank[li3|D]={rankM}")
def wedge(u,v):
    out={}
    for a,ca in u.items():
        for b,cb in v.items():
            if a==b: continue
            key,sg=((a,b),1) if a<b else ((b,a),-1)
            out[key]=out.get(key,Fr(0))+sg*ca*cb
    return {k:c for k,c in out.items() if c}
def delta_li3(z):
    vz=S.vec(z); vomz=S.vec(E.fsub(E.ONE,z)); w2=wedge(vz,vomz); out={}
    for ab,cab in w2.items():
        for l3,c3v in vz.items(): out[(ab,l3)]=out.get((ab,l3),Fr(0))+cab*c3v
    return {k:c for k,c in out.items() if c}
DL=[delta_li3(z) for z in CRs]
Mnull=M.nullspace(); Rli3=[]
for vec in Mnull:
    cc=[Fr(sp.Rational(vec[z])) for z in range(24)]
    if any(x!=0 for x in cc): Rli3.append(cc)
def delta_of(c):
    acc={}
    for z,cz in enumerate(c):
        if cz:
            for key,cf in DL[z].items(): acc[key]=acc.get(key,Fr(0))+cf*cz
    return {k:v for k,v in acc.items() if v}
rel_ok=all(not delta_of(r) for r in Rli3)
check("B4: delta vanishes on all li3-relations (well-defined)", rel_ok, f"{len(Rli3)} relations")
dcomp=[{},{},{}]
for k in range(3):
    free=set()
    for z in range(sols[k].rows): free|=sols[k][z].free_symbols
    ck=[Fr(sp.Rational(sols[k][z].subs({fs:0 for fs in free}))) for z in range(24)]
    for z in range(24):
        if ck[z]:
            for key,cf in DL[z].items(): dcomp[k][key]=dcomp[k].get(key,Fr(0))+cf*ck[z]
allk=set(dcomp[0])|set(dcomp[1])|set(dcomp[2])
delta={}
for key in allk:
    A=dcomp[0].get(key,Fr(0)); B=dcomp[1].get(key,Fr(0)); C=dcomp[2].get(key,Fr(0))
    if A or B or C: delta[key]=(A,B,C,1)
nonzero = len(delta)>0
exp = {((1,2),1):(Fr(601,32),Fr(-91,16),Fr(-243,32)), ((1,2),2):(Fr(613,32),Fr(-59,16),Fr(-73,8))}
match = all(delta.get(k,(0,0,0,1))[:3]==v for k,v in exp.items())
check("B5: delta(N3) NONZERO", nonzero, f"{len(delta)} terms")
check("B6: delta matches reported exact coefficients", match, str({k:v[:3] for k,v in sorted(delta.items())}))
print("\nOVERALL:", "PASS" if PASS[0] else "FAIL")
'''
open('verify_pg6.py','w').write(verify)
print("verify_pg6.py written OK,", len(verify), "chars")
