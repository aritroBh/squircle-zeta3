import sys, json, math, time, itertools
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import flint, mpmath as mp
import w5blib as W, w5verify as V

gl = json.load(open('out/d13/genuine_lattice_cap.json'))
G = [[int(x) for x in row] for row in gl['G']]
atom_names = gl['atom_names']; n=len(atom_names)
J = json.load(open('out/d13/candidates_cap.json'))
Kk = json.load(open('out/d13/kernel5_cap.json'))
K=Kk['K']; comps=Kk['comps']; names=J['names']
keepc = [0,1,2,3]+list(range(23,51))
g = 36771482407
def canon_key(parts): return tuple(sorted(parts))
prodkeys=set()
for kj in keepc:
    for ms_s in comps[kj]:
        prodkeys.add(canon_key(names[ii] for ii in eval(ms_s)))
for pname in atom_names[2+len(keepc):]:
    prodkeys.add(canon_key(tuple(pname.split('*'))))
prodkeys=sorted(prodkeys); PK={k:i for i,k in enumerate(prodkeys)}
L=120
def expand(v): return V.expand_a(13, J['cands'], K, comps, keepc, atom_names, v, names)
def omap(v):
    e=expand(v)
    li5=[0]*len(J['cands'])
    for i,c in e['li5'].items(): li5[i]=int(c)*L
    pr=[0]*len(prodkeys)
    for k,c in e['prod'].items(): pr[PK[canon_key(k)]]+=int(c*L)
    return li5+pr+[int(e['a_z5'])*L]
Ts=[row[0] for row in G]
ker = W.int_left_kernel_hnf([[t] for t in Ts])
L0=[]
for cvec in ker:
    r=[0]*n
    for c,i in zip(cvec, range(9)):
        if c:
            for j in range(n): r[j]+=c*G[i][j]
    L0.append(r)
R=flint.fmpz_mat(L0).lll(gram='exact')
L0r=[[int(R[i,j]) for j in range(n)] for i in range(R.nrows()) if any(int(R[i,j]) for j in range(n))]
vatom=[int(x) for x in json.load(open('out/d13/gcd_vec_atom.json'))['v_atom']]
Omat=[omap(r) for r in L0r]; ncols=len(Omat[0])
Ro=flint.fmpz_mat(Omat).lll(gram='exact')
Ror=[[int(Ro[i,j]) for j in range(ncols)] for i in range(Ro.nrows())]
import sympy as sp
Om=sp.Matrix(Omat)
piv=list(Om.rref()[1])
print('Om rank pivots:', piv)
assert len(piv)==8
C=sp.Matrix([[Omat[i][j] for j in piv] for i in range(8)])
RC=sp.Matrix([[Ror[i][j] for j in piv] for i in range(8)])
Ci=C.inv()
U=RC*Ci
assert all(x.q==1 for x in U), 'U not integral'
U=[[int(x) for x in row] for row in U.tolist()]
detU=sp.Matrix(U).det()
print('U det:', detU)
assert abs(detU)==1
# verify U*Omat == Ror
chk=(sp.Matrix(U)*Om - sp.Matrix(Ror)).norm()
print('transform check (should be 0):', chk)
L0o=[[sum(U[i][kk]*L0r[kk][j] for kk in range(8)) for j in range(n)] for i in range(8)]
for r in L0o: assert r[0]==0
json.dump(dict(L0o=[[str(x) for x in r] for r in L0o], L=str(L), prodkeys=['*'.join(k) for k in prodkeys]),
          open('out/d13/gcd_L0_ometric.json','w'))
print('saved; L0o atom-max per row:', [max(abs(x) for x in r) for r in L0o])

# Babai of O(vatom) against reduced O-basis
mp.mp.dps=300
B=[[mp.mpf(x) for x in r] for r in Ror]
m=ncols; kk=8
Bs=[]; mu=[[mp.mpf(0)]*kk for _ in range(kk)]; nrm=[]
for i in range(kk):
    bi=B[i][:]
    for j in range(i):
        mu[i][j]=sum(B[i][t]*Bs[j][t] for t in range(m))/sum(Bs[j][t]**2 for t in range(m))
        for t in range(m): bi[t]-=mu[i][j]*Bs[j][t]
    Bs.append(bi); nrm.append(sum(x*x for x in bi))
Ov=omap(vatom)
b=[mp.mpf(x) for x in Ov]; coeff=[0]*kk
for i in range(kk-1,-1,-1):
    c=sum(b[t]*Bs[i][t] for t in range(m))/nrm[i]
    ci=int(mp.nint(c)); coeff[i]=ci
    fr=abs(float(c-ci))
    if fr>0.49: print('  WARN half-int at', i, fr)
    for t in range(m): b[t]-=ci*B[i][t]
print('Babai O-coeffs:', coeff)
vnum=[vatom[j]-sum(coeff[i]*L0o[i][j] for i in range(kk)) for j in range(n)]
assert vnum[0]==g, vnum[0]
e=expand(vnum)
def stats(e):
    coefs=list(e['li5'].values())+list(e['prod'].values())+[e['a_z5']]
    mn=max(abs(c.numerator) for c in coefs); lm=1
    for c in coefs: lm=math.lcm(lm,c.denominator)
    return dict(nli5=len(e['li5']), nprod=len(e['prod']), maxnum=mn, lcm=lm, z5=str(e['a_z5']))
print('v_num (Babai-O):', stats(e), 'atom-max:', max(abs(x) for x in vnum))
json.dump(dict(v_num=[str(x) for x in vnum], babai_coeff=coeff), open('out/d13/gcd_vec_num_babai.json','w'))
