import os, sys, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
# check flint nmod_mat availability
import flint
print(flint.__version__ if hasattr(flint,'__version__') else 'flint')
try:
    from flint import nmod_mat
    M=nmod_mat(3,3,[1,2,3,4,5,6,7,8,9],17)
    print("nmod_mat rank:", M.rank())
    HAS_NMOD=True
except Exception as e:
    print("nmod_mat:", e); HAS_NMOD=False
# build degree-12 monomials + reversal orbits
def deg_monomials(d, n=5):
    out=[]
    def rec(i, rem, cur):
        if i==n-1: out.append(tuple(cur+[rem])); return
        for k in range(rem+1): rec(i+1, rem-k, cur+[k])
    rec(0,d,[]); return out
mons12=deg_monomials(12)
print("deg-12 monomials:", len(mons12), "(expect C(16,4)=1820)")
seen=set(); orbits12=[]
for e in mons12:
    if e in seen: continue
    r=e[::-1]; orbits12.append((e,) if r==e else (e,r)); seen.add(e); seen.add(r)
print("reversal orbits deg-12:", len(orbits12), "(expect 924)")
# degree-8 for quotients
mons8=deg_monomials(8)
seen=set(); orbits8=[]
for e in mons8:
    if e in seen: continue
    r=e[::-1]; orbits8.append((e,) if r==e else (e,r)); seen.add(e); seen.add(r)
print("deg-8 monomials:", len(mons8), "orbits:", len(orbits8), "(expect 495, 255)")
