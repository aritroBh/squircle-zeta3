import os, json, time, sys
cwd='<project>/research/08_theory/DCYA1_corner_decay_proof/'
os.chdir(cwd)
sys.path.insert(0,cwd)
print("Files present:", [f for f in sorted(os.listdir(cwd)) if f.endswith('.json') or f.endswith('.py') or f.endswith('.md')])
import mpmath as mp
mp.mp.dps=30
# Check tables
for f in ['mk_2pi3_k80.json','mk_pi2_k100_dps140.json']:
    if os.path.exists(f):
        d=json.load(open(f)); print(f, "has", len(d), "entries; first:", d[0][:20])
# Confirm alpha=2pi/3 constants
a=2*mp.pi/3
A=mp.sqrt(2*a*mp.sin(a))
print("alpha=2pi/3: A=",mp.nstr(A,10)," sqrt(2a)=",mp.nstr(mp.sqrt(2*a),8)," delta_rest=sqrt(2(2pi-a))=",mp.nstr(mp.sqrt(2*(2*mp.pi-a)),8))
print("delta_rest > sqrt(2a)?", mp.sqrt(2*(2*mp.pi-a))>mp.sqrt(2*a))
