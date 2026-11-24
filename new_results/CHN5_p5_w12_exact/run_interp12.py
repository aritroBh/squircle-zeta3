
import sys, time, pickle
sys.path.insert(0,'.')
from fractions import Fraction as F
from linalg_tools import solve_linear
import sympy as sp
dd=pickle.load(open('det12_points.pkl','rb'))
sel=[tuple(x) for x in dd['sel']]
orbits12=dd['orbits12']; mons12=dd['mons12']; orbits8=dd['orbits8']; mons8=dd['mons8']
interp_idx=dd['interp_idx']
# load all eval chunks
import glob, re
chunks={}
for fn in glob.glob('eval12_*.pkl'):
    st=int(re.search(r'eval12_(\d+)\.pkl',fn).group(1))
    chunks[st]=pickle.load(open(fn,'rb'))
data=[]
for st in sorted(chunks):
    data.extend(chunks[st])
assert len(data)==len(sel), (len(data),len(sel))
assert [tuple(d['pt']) for d in data]==[tuple(x) for x in sel], "point order mismatch"
N8=len(orbits8); N12=len(orbits12)
orbit_of8={}
for i,o in enumerate(orbits8):
    for e in o: orbit_of8[e]=i
def orbit_vals(a, mons, orbit_of, N):
    vals=[F(0)]*N
    for e in mons:
        pw=F(1)
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        vals[orbit_of[e]]+=pw
    return vals
def Lval(a):
    a1,a2,a3,a4,a5=a; return F(a2)*a3*a4*(a1+a2+a3+a4+a5)
# coordinate extractors
def coords_of(d):
    g=d['g']; m0=d['m0']
    ft=[g[i]+m0[i] for i in range(12)]
    return {'A':ft[7],'B':ft[8],'CLASS':ft[9]+F(14,9)*ft[10],'Z':ft[11],
            'm0_A':m0[7],'m0_B':m0[8],'m0_CLASS':m0[9]+F(14,9)*m0[10],'m0_Z':m0[11]}
rows8={}
with open('interp12.log','w') as lg:
    # 1) check m0 vanishing at all 924
    maxm0={k:F(0) for k in ['m0_A','m0_B','m0_CLASS','m0_Z']}
    for d in data:
        c=coords_of(d)
        for k in maxm0: maxm0[k]=max(maxm0[k],abs(c[k]))
    lg.write(f"m0 max |contribution| over 924 pts: {maxm0}\n"); lg.flush()
    # 2) interpolate each quotient
    results={}
    for X in ['A','B','CLASS','Z']:
        rows=[orbit_vals(data[i]['pt'],mons8,orbit_of8,N8) for i in interp_idx]
        b=[coords_of(data[i])[X]/Lval(data[i]['pt']) for i in interp_idx]
        t0=time.time()
        sol=solve_linear(rows,b)
        if sol is None:
            lg.write(f"{X}: INCONSISTENT interpolation system!\n"); lg.flush(); continue
        Q={}
        for oi,cf in zip(range(N8),sol):
            if cf:
                for e in orbits8[oi]:
                    Q[e]=Q.get(e,F(0))+cf
        lg.write(f"{X}: interpolated quotient, {sum(1 for v in sol if v)} nonzero orbit coeffs, {time.time()-t0:.0f}s\n"); lg.flush()
        results[X]=Q
    pickle.dump(results, open('Q8_quotients.pkl','wb'))
    # 3) verify D = f - L*Q8 == 0 at all 924 points
    def peval_orbit(Q, a):
        tot=F(0)
        for e,c in Q.items():
            pw=F(1)
            for x,ei in zip(a,e):
                if ei: pw*=x**ei
            tot+=c*pw
        return tot
    for X in ['A','B','CLASS','Z']:
        if X not in results: continue
        Q=results[X]; bad=0; t0=time.time()
        for d in data:
            D=coords_of(d)[X]-Lval(d['pt'])*peval_orbit(Q,d['pt'])
            if D!=0: bad+=1; lg.write(f"  {X} NONZERO D at {d['pt']}: {D}\n")
        lg.write(f"{X}: D==0 at all {len(data)} pts: {bad==0} (bad={bad}) [{time.time()-t0:.0f}s]\n"); lg.flush()
    lg.write("DONE\n")
print("interp done")
