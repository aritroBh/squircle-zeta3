import os
cwd = os.getcwd()
gscript = r'''import sys, json, time, pickle, math
from fractions import Fraction as Ff
import mpmath as mp
import flint
sys.path.insert(0,'.')
import w5blib as W
raw = pickle.load(open('out/d13/atoms5_cap_d17.pkl','rb'))
atom_names=raw['atom_names']
mp.mp.dps=2050
atom_vals=[mp.mpf(v) for v in raw['atom_vals']]
n=len(atom_vals)
P=1800
mp.mp.dps=P+30
Cc=mp.power(10,P)
B=[]
for i,v in enumerate(atom_vals):
    row=[0]*(n+1); row[i]=1; row[n]=int(mp.nint(v*Cc)); B.append(row)
t1=time.time(); R=flint.fmpz_mat(B).lll(gram='exact'); dt=time.time()-t1
print(f"LLL done {dt:.0f}s", flush=True)
mp.mp.dps=2050
genuine=[]
for i in range(n):
    a=[int(R[i,j]) for j in range(n)]
    if not any(a): continue
    res=sum(mp.mpf(x)*v for x,v in zip(a,atom_vals))
    genuine.append((abs(res), a, res))
genuine.sort(key=lambda t: t[0])
# classify by residual
print("residual profile (smallest 8):", [mp.nstr(g[0],2) for g in genuine[:8]], flush=True)
print("residual profile (around rank):", [mp.nstr(g[0],2) for g in genuine[60:74]], flush=True)
# genuine = residual < 1e-1900 (atom floor ~1e-2008, spurious ~1e-1800)
G=[a for r,a,rs in genuine if r < mp.mpf('1e-1900')]
print(f"genuine relations (resid<1e-1900): {len(G)} of {len(genuine)}", flush=True)
# T-coef gcd over genuine lattice
g=0
for a in G: g=math.gcd(g, abs(a[0]))
print(f"gcd of T-coefs over genuine lattice: {g}", flush=True)
# min |a[0]| nonzero
mn=min((abs(a[0]) for a in G if a[0]), default=0)
print(f"min |T-coef| among LLL rows: {mn}", flush=True)
# HNF on T-coef column to find min achievable and vector
import sympy as sp
M=sp.Matrix(G)
# Compute integer relation for min T-coef via HNF on augmented [a0 | rest]
# Save genuine lattice for minimization
json.dump(dict(G=[[str(x) for x in a] for a in G], atom_names=atom_names),
          open('out/d13/genuine_lattice_cap.json','w'))
print("saved out/d13/genuine_lattice_cap.json", flush=True)
# T-free sublattice rank
G0=[a for a in G if a[0]==0]
print(f"T-free genuine relations: {len(G0)}", flush=True)
'''
open('step_24_d13_lattice.py','w').write(gscript)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_24_d13_lattice.py'],
                   capture_output=True, text=True, cwd=cwd, timeout=2400)
print("elapsed", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[-1500:] if r.stderr else "")
