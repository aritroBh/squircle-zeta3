proof = r'''
import sys, time, pickle
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F

def make_f_eval(W, idx):
    T = PE.table(W)
    def f(a):
        a=[F(x) for x in a]
        acc=F(0)
        for pt in (a, a[::-1]):
            vw = PE.volreg_words_point(pt, W)
            for w,pp in vw[W].items():
                c=T.reduce_word(w)[idx]
                if c: acc+=pp*c
        return acc
    return f
f11 = make_f_eval(11,8)
_ = f11((1,2,3,4,5))

dd=pickle.load(open('det11_points.pkl','rb')); sel=[tuple(x) for x in dd['sel']; orbits11=dd['orbits11']; mons11=dd['mons11']
Q11=pickle.load(open('Q11_poly.pkl','rb'))
def Lval(a):
    a1,a2,a3,a4,a5=a; return a2*a3*a4*(a1+a2+a3+a4+a5)
def peval(poly,a):
    a=[F(x) for x in a]; tot=F(0)
    for e,c in poly.items():
        pw=F(1)
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        tot+=c*pw
    return tot
# rank check of det11 set (deg-11 reversal orbits) mod 2 primes
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
rows=[orbit_vals11(a) for a in sel]
with open('proof_div.log','w') as lg:
    for p in [(1<<61)-1, (1<<63)-25]:
        rk=rank_mod(rows,p)
        lg.write(f"det11 rank mod {p}: {rk} (need {N11})\n"); lg.flush()
    assert all(rank_mod(rows,p)==N11 for p in [(1<<61)-1]), "rank not full"
    # evaluate D = f11 - L*Q11 at each point
    bad=0; t0=time.time()
    for i,a in enumerate(sel):
        D=f11(a)-Lval(a)*peval(Q11,a)
        if D!=0: bad+=1; lg.write(f"NONZERO at {a}: {D}\n")
        if (i+1)%25==0: lg.write(f"{i+1}/{len(sel)} checked, bad={bad}, {time.time()-t0:.0f}s\n"); lg.flush()
    lg.write(f"DONE: {len(sel)} points, {bad} nonzeros, {time.time()-t0:.0f}s\n")
    lg.write("VERDICT: "+("f11 = a2a3a4(Sum a)Q11 PROVED (D==0 on rank-693 determining set)" if bad==0 else "FAIL")+"\n")
print("proof done")
'''
open('run_proof_div.py','w').write(proof)
import subprocess
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_proof_div.py'],
                      stdout=open('run_proof_div.out','w'), stderr=subprocess.STDOUT)
open('run_proof_div.pid','w').write(str(proc.pid))
print("launched proof pid", proc.pid)
