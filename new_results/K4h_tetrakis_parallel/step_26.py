step04b = r'''"""step_04b_guard.py -- decompose the val96 diff: parallel at dps 80 with guard=50
(workers accumulate at 130 digits) -> S_g ~ the exact-arithmetic value of the
same quadrature rule. Then E_par = |S_g - S_par| (parallel scheme error) and
E_ref = |S_g - S_ref| (single-process chamber_S error) at dps 80.
"""
import sys, os, time, json
sys.dont_write_bytecode = True
os.chdir("<project>/research/08_theory/K4h_tetrakis_parallel")
K4E = "<project>/research/08_theory/K4e_tetrakis_d3"
if K4E not in sys.path: sys.path.insert(0, K4E)
import mpmath as mp
from mpmath import mpf

LOGF = "val96g.log"
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    open(LOGF, "a").write(s + "\n")

def main():
    mp.mp.dps = 80
    import k4h_driver as KD
    s2 = mp.sqrt(2)
    NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
    H = mpf(1)/96; T = mpf('4.4')
    cfg = dict(normals=[[str(x) for x in row] for row in NA], K=4,
               HO=str(H), HI=str(H), HI2=str(H), TO=str(T), TI=str(T),
               WMIN='1e-45', LMIN='1e-35', dps=80, guard=50, chunk=32)
    open(LOGF,"w").write("val96 guarded parallel (workers at dps 130)\n")
    t0 = time.time()
    acc, stats, N, nch = KD.run_job(cfg, "ckpt_val96_g130", 6, log, "val96g")
    S_g = [48*v for v in acc]
    log("[val96g] guarded parallel done: %.0fs N=%d nch=%d stats=%s" % (time.time()-t0, N, nch, stats))

    prev = json.load(open("val96_result.json"))
    mp.mp.dps = 140
    S_par = [mpf(s) for s in prev["Spar"]]
    S_ref = [mpf(s) for s in prev["Sref"]]
    print()
    for k in range(5):
        ep = abs(S_g[k] - S_par[k]); er = abs(S_g[k] - S_ref[k]); dpr = abs(S_par[k]-S_ref[k])
        log("k=%d  E_par=|Sg-Spar|=%s   E_ref=|Sg-Sref|=%s   |par-ref|=%s"
            % (k, mp.nstr(ep,3), mp.nstr(er,3), mp.nstr(dpr,3)))
    log("[val96g] DONE")

if __name__ == "__main__":
    main()
'''
open("step_04b_guard.py","w").write(step04b)
import subprocess, sys, os, time
subprocess.run(["rm","-rf","ckpt_val96_g130"],capture_output=True)
out = open("val96g_out.txt","w")
p = subprocess.Popen(["nohup", sys.executable, "step_04b_guard.py"], stdout=out, stderr=subprocess.STDOUT,
                     stdin=subprocess.DEVNULL, start_new_session=True)
print("step_04b launched, PID", p.pid)
time.sleep(30)
print(open("val96g_out.txt").read())
