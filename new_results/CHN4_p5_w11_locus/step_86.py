import time, pickle
def log_progress(line):
    with open('progress.log','a') as f: f.write(time.strftime('%H:%M:%S ')+line+'\n')
log_progress("DIVISIBILITY PROVED: D=f11-a2a3a4(Sum a)Q11 = 0 at all 693 points of rank-693 reversal-invariant deg-11 determining set (run_proof_div.py, 1574s, 0 nonzeros). => f11 = -(1/10)a2a3a4(Sum a)P7 EXACTLY. Upgrades Stage2 interpolation from 48-pt verification to full PROOF.")
# re-verify det11 rank mod 2 primes (proof depends on it)
dd=pickle.load(open('det11_points.pkl','rb')); sel=dd['sel']; orbits11=dd['orbits11']; mons11=dd['mons11']
orbit_of11={}
for i,o in enumerate(orbits11):
    for e in o: orbit_of11[e]=i
N11=len(orbits11)
def orbit_vals11(a):
    vals=[0]*N11
    for e in mons11:
        pw=1
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        vals[orbit_of11[e]]+=pw
    return vals
def rank_mod(rows,p):
    M=[[x%p for x in r] for r in rows]; nr=len(M); nc=len(M[0]); r=0
    for c in range(nc):
        piv=next((i for i in range(r,nr) if M[i][c]%p),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c]%p,p-2,p); M[r]=[(v*inv)%p for v in M[r]]
        for i in range(nr):
            if i!=r and M[i][c]%p:
                f=M[i][c]%p; M[i]=[(M[i][j]-f*M[r][j])%p for j in range(nc)]
        r+=1
        if r==nr: break
    return r
import time
t0=time.time()
rows=[orbit_vals11(a) for a in sel]
for p in [(1<<61)-1, (1<<63)-25]:
    print(f"det11 rank mod {p}: {rank_mod(rows,p)} (need {N11})  [{time.time()-t0:.0f}s]")
