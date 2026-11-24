import sys, subprocess
lib = open('w5lib.py').read()
lib = lib.replace("""    # sum over k>=0, k!=4: zeta(5-k)(-t)^k/k!
    k=0
    term_thresh=mp.mpf(10)**(-(P+15))
    while k<300:
        if k!=4:
            zk=mp.zeta(5-k)
            term=zk*(-t)**k/mp.factorial(k)
            s+=term
            if k>12 and abs(term)<term_thresh:
                break
        k+=1""",
"""    # sum over k>=0, k!=4: zeta(5-k)(-t)^k/k!   (fixed count; terms with zeta(even neg)=0 are skipped)
    for k in range(0, 320):
        if k==4: continue
        zk=mp.zeta(5-k)
        if zk:
            s+=zk*(-t)**k/mp.factorial(k)""")
open('w5lib.py','w').write(lib)
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
print("jonquiere check v2:", r0.stdout, r0.stderr[-300:])
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_49_d21_lll.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_49_d21_lll(r3): jonquiere fix; D=21 numerics+LLL rerun\n")
