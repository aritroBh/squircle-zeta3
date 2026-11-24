import sys, json, math, time, itertools
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import flint
import w5blib as W, w5verify as V

SRC = '''# step_01_gcd.py -- see transcript step file (code identical): gcd combination + T-free size reduction'''
d='.'
G = [[int(x) for x in row] for row in json.load(open('out/d13/genuine_lattice_cap.json'))['G']]
gl = json.load(open('out/d13/genuine_lattice_cap.json'))
atom_names = gl['atom_names']
J = json.load(open('out/d13/candidates_cap.json'))
Kk = json.load(open('out/d13/kernel5_cap.json'))
names = J['names']
keepc = [0,1,2,3]+list(range(23,51))
assert atom_names == ['T13^5','z5']+[f'K{j}' for j in keepc]+atom_names[2+len(keepc):]
n=len(atom_names)

# --- (a) T column, gcd, extended gcd
Ts = [row[0] for row in G]
g = 0
for t in Ts: g = math.gcd(g, t)
print('T coefs:', Ts)
print('gcd =', g, ' factor:', __import__('sympy').factorint(g))
assert g == 36771482407
# extended gcd over the 5 nonzero T rows (idx 4..8)
def xgcd_list(vals):
    g0=abs(vals[0]); cs=[1 if vals[0]>=0 else -1]
    for v in vals[1:]:
        # combine g0 and v
        a,b=g0,v
        s0,s1,t0,t1=1,0,0,1
        while b:
            q=a//b; a,b=b,a-q*b; s0,s1=s1,s0-q*s1; t0,t1=t1,t0-q*t1
        # a = s0*g0 + t0*v ; keep sign positive
        if a<0: a=-a; s0=-s0; t0=-t0
        cs=[c*s0 for c in cs]+[t0]; g0=a
    return g0, cs
idxT=[i for i,t in enumerate(Ts) if t]
g2, cs = xgcd_list([Ts[i] for i in idxT])
assert g2==g
v0=[0]*n
for c,i in zip(cs, idxT):
    for j in range(n): v0[j]+=c*G[i][j]
assert v0[0]==g
print('Bezout coeffs:', dict(zip(idxT,cs)), ' v0 max|coef|:', max(abs(x) for x in v0))

# --- (b) T-free sublattice: left kernel of T column (rank 8)
M=[[t] for t in Ts]  # 9x1
ker = W.int_left_kernel_hnf(M)   # list of 9-vectors, c.Ts=0
print('T-free kernel rank:', len(ker))
L0=[]
for cvec in ker:
    r=[0]*n
    for c,i in zip(cvec, range(9)):
        if c:
            row=G[i]
            for j in range(n): r[j]+=c*row[j]
    L0.append(r)
t0=time.time()
R=flint.fmpz_mat(L0).lll(gram='exact')
L0r=[[int(R[i,j]) for j in range(n)] for i in range(R.nrows())]
L0r=[r for r in L0r if any(r)]
print(f'LLL on T-free lattice {len(L0r)}x{n} done ({time.time()-t0:.0f}s); row max coefs:', [max(abs(x) for x in r) for r in L0r])

# --- (c) Babai nearest-plane (mpmath 200 digits) then enumerate {-1,0,1}^8
import mpmath as mp
mp.mp.dps=200
B=[[mp.mpf(x) for x in r] for r in L0r]
k=len(B)
# GSO
Bs=[]; mu=[[mp.mpf(0)]*k for _ in range(k)]
for i in range(k):
    bi=B[i][:]
    for j in range(i):
        mu[i][j]=sum(B[i][t]*Bs[j][t] for t in range(n))/sum(Bs[j][t]**2 for t in range(n))
        for t in range(n): bi[t]-=mu[i][j]*Bs[j][t]
    Bs.append(bi)
def babai(v):
    b=[mp.mpf(x) for x in v]; coeff=[0]*k
    for i in range(k-1,-1,-1):
        c=sum(b[t]*Bs[i][t] for t in range(n))/sum(Bs[i][t]**2 for t in range(n))
        ci=int(mp.nint(c)); coeff[i]=ci
        frac=abs(c-ci)
        if frac>mp.mpf('0.49') and frac<mp.mpf('0.51'): print('  WARN near half-integer:', i, mp.nstr(frac,5))
        for t in range(n): b[t]-=ci*B[i][t]
    return coeff, b
coeff, rem = babai(v0)
v1=[v0[j]-sum(coeff[i]*L0r[i][j] for i in range(k)) for j in range(n)]
assert v1[0]==g
print('Babai: shift coeffs', coeff, ' v1 max|coef|:', max(abs(x) for x in v1))
# enumerate around Babai
best=None
cnt=0
t0=time.time()
for shift in itertools.product((-1,0,1), repeat=k):
    vv=[v1[j]+sum(shift[i]*L0r[i][j] for i in range(k)) for j in range(n)]
    mx=max(abs(x) for x in vv)
    l2=sum(x*x for x in vv)
    cnt+=1
    if best is None or (mx,l2)<(best[0],best[1]): best=(mx,l2,vv[:])
print(f'enum {cnt} shifts ({time.time()-t0:.0f}s); best atom max|coef| {best[0]}, L2 {best[1]}')
vatom=best[2]
assert vatom[0]==g
json.dump(dict(v_atom=[str(x) for x in vatom], atom_max=best[0]),
          open('out/d13/gcd_vec_atom.json','w'))
print('saved out/d13/gcd_vec_atom.json')
print('nonzero entries of v_atom:', sum(1 for x in vatom if x))
