import sys, subprocess
lib = open('w5lib.py').read()
lib = lib.replace("""def li5_series(x,P):""", """def li5_jonquiere(x,P):
    \"\"\"Li5(x) for 0<x<1 near 1 via Jonquiere's expansion, t=-log(x), |t|<2pi.\"\"\"
    old=mp.mp.dps; mp.mp.dps=P+30
    x=mp.mpf(x)
    t=-mp.log(x)
    H4=mp.mpf(25)/12
    s=t**4/24*(H4-mp.log(t))
    # sum over k>=0, k!=4: zeta(5-k)(-t)^k/k!
    k=0
    term_thresh=mp.mpf(10)**(-(P+15))
    while k<300:
        if k!=4:
            zk=mp.zeta(5-k)
            term=zk*(-t)**k/mp.factorial(k)
            s+=term
            if k>12 and abs(term)<term_thresh:
                break
        k+=1
    mp.mp.dps=old
    return s

def li5_smart(x,P):
    x=mp.mpf(x)
    if mp.mpf('0.9') < x < 1:
        return li5_jonquiere(x,P)
    return li5_series(x,P)

def li5_series(x,P):""")
open('w5lib.py','w').write(lib)
p = open('pipe.py').read()
p = p.replace("        v = W.li5_series(x, P)", "        v = W.li5_smart(x, P)")
open('pipe.py','w').write(p)
# cross-check jonquiere vs series at several points
chk = r'''
import sys; sys.path.insert(0,'.')
import mpmath as mp
import w5lib as W
for xs in (0.95, 0.9909, 0.999917, 0.970562748477141):
    mp.mp.dps = 300
    v1 = W.li5_smart(mp.mpf(xs), 250)
    v2 = mp.polylog(5, mp.mpf(xs))
    print(xs, mp.nstr(abs(v1-v2), 3))
'''
r0 = subprocess.run(['nice','-n','10',sys.executable,'-c',chk],capture_output=True,text=True,timeout=900)
print("jonquiere check:", r0.stdout, r0.stderr[-300:])
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_49_d21_lll.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_49_d21_lll(r2): li5_smart; D=21 numerics + kernel LLL\n")
