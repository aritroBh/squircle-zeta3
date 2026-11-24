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
# product key universe (canonical)
prodkeys=set()
for kj in keepc:
    for ms_s in comps[kj]:
        ms=eval(ms_s); prodkeys.add(canon_key(names[ii] for ii in ms))
for pname in atom_names[2+len(keepc):]:
    prodkeys.add(canon_key(tuple(pname.split('*'))))
prodkeys=sorted(prodkeys)
PK={k:i for i,k in enumerate(prodkeys)}
# L = lcm of all companion denominators in keepc
L=1
for kj in keepc:
    for cf_s in comps[kj].values(): L=math.lcm(L, Ff(cf_s).denominator)
print('companion denom lcm L =', L, '; product keys:', len(prodkeys))

def expand(v):
    return V.expand_a(13, J['cands'], K, comps, keepc, atom_names, v, names)

def omap(v):
    e=expand(v)
    li5=[0]*len(J['cands'])
    for i,c in e['li5'].items(): li5[i]=int(c)*L
    pr=[0]*len(prodkeys)
    for k,c in e['prod'].items():
        ck=canon_key(k); pr[PK[ck]]+= int(c*L)
    return li5+pr+[int(e['a_z5'])*L]

# T-free lattice L0r (atom space) from previous step: rebuild
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
assert vatom[0]==g

# baseline expansion stats of v_atom
e0=expand(vatom)
def stats(e):
    coefs=list(e['li5'].values())+list(e['prod'].values())+[e['a_z5']]
    mn=max(abs(c.numerator) for c in coefs)
    lm=1
    for c in coefs: lm=math.lcm(lm,c.denominator)
    return dict(nli5=len(e['li5']), nprod=len(e['prod']), maxnum=mn, lcm=lm, z5=str(e['a_z5']))
print('v_atom expansion:', stats(e0))

# build O-matrix of T-free basis (8 x (275+P+1))
t0=time.time()
Omat=[omap(r) for r in L0r]
ncols=len(Omat[0])
print(f'Omat built {len(Omat)}x{ncols} ({time.time()-t0:.0f}s); max entry {max(abs(x) for r in Omat for x in r)}')
t0=time.time()
Ro=flint.fmpz_mat(Omat).lll(gram='exact')
Ror=[[int(Ro[i,j]) for j in range(ncols)] for i in range(Ro.nrows())]
print(f'O-metric LLL done ({time.time()-t0:.0f}s); row norms:', [int(mp.sqrt(sum(mp.mpf(x)**2 for x in r))) for r in Ror])
# recover transform U = R * Omat^+ using 8 independent columns
import sympy as sp
Om=sp.Matrix(Omat)
piv = list(Om.T.rref()[1])[:8]  # pivot columns? rref of transpose gives pivot rows = indep columns
print('pivot cols:', piv)
C=sp.Matrix([[Omat[i][j] for j in piv] for i in range(8)])
RC=sp.Matrix([[Ror[i][j] for j in piv] for i in range(8)])
U=RC*C.inv()
assert all(x.q==1 for x in U), 'U not integral?!'
U=[[int(x) for x in row] for row in U.tolist()]
print('U det:', abs(round(sp.Matrix(U).det())), 'U:', U)
# new atom-space T-free basis with reduced O-image
L0o=[[sum(U[i][kk]*L0r[kk][j] for kk in range(8)) for j in range(n)] for i in range(8)]
json.dump(dict(L0o=[[str(x) for x in r] for r in L0o], L=str(L), prodkeys=['*'.join(k) for k in prodkeys]),
          open('out/d13/gcd_L0_ometric.json','w'))
print('saved out/d13/gcd_L0_ometric.json; L0o row atom-max:', [max(abs(x) for x in r) for r in L0o])
