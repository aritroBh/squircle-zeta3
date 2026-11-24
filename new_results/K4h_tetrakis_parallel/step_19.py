job = r'''#!/usr/bin/env python3
"""jobA_h768.py -- K4h production: model A (tetrakis hexahedron), h=1/768,
dps 270, WMIN 1e-215, LMIN 1e-205, TO=TI=6.0, 6 workers, chunk=24,
checkpoints in ckpt_h768/ (resume: just rerun this script).

Parallel equivalent of K4e's jobA_h384.py: identical quadrature, split over
outer phi-nodes (k4h_driver.run_job). Writes modelA_h768.json in the same
format as modelA_h384.json: {"Sdk": [str x5], "stats": {...}, "HO": str, "HI": str}.
Provenance extras go to modelA_h768_meta.json.
"""
import os, sys, json, time
CWD = "<project>/research/08_theory/K4h_tetrakis_parallel"
os.chdir(CWD)
sys.dont_write_bytecode = True
K4E = "<project>/research/08_theory/K4e_tetrakis_d3"
if K4E not in sys.path: sys.path.insert(0, K4E)
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 270
import k4h_driver as KD

DPS = 270; K = 4
H = mpf(1)/768
T = mpf('6.0')
WMIN = mpf('1e-215'); LMIN = mpf('1e-205')
WORKERS = 6; CHUNK = 24
CKDIR = "ckpt_h768"; LOG = "bgA768.log"

def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    open(LOG, "a").write(s + "\n")

def main():
    s2 = mp.sqrt(2)
    NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
    open("jobA_h768.pid", "w").write(str(os.getpid()))
    open(LOG, "a").write("jobA768 start dps=%d pid=%d nice=%d time=%s\n"
                         % (DPS, os.getpid(), os.nice(0), time.strftime("%Y-%m-%d %H:%M:%S")))
    cfg = dict(normals=[[str(x) for x in row] for row in NA], K=K,
               HO=str(H), HI=str(H), HI2=str(H), TO=str(T), TI=str(T),
               WMIN=str(WMIN), LMIN=str(LMIN), dps=DPS, chunk=CHUNK)
    t0 = time.time()
    acc, stats, N, nch = KD.run_job(cfg, CKDIR, WORKERS, log, "h768")
    el = time.time() - t0
    S = [48*v for v in acc]
    S1 = 2*mp.sqrt(5)/3
    vals = []
    for k in range(K+1):
        v = mp.re(S[k])
        vals.append(mp.nstr(v, DPS))
        log("[h768] [d^%d]S = %s" % (k, mp.nstr(v, 250)))
    log("[h768] el=%ds stats=%s" % (el, stats))
    log("[h768] d0 err vs 2sqrt5/3 = %s" % mp.nstr(abs(mp.re(S[0])-S1), 3))
    json.dump({"Sdk": vals, "stats": stats, "HO": str(H), "HI": str(H)},
              open("modelA_h768.json", "w"), indent=1)
    json.dump({"dps": DPS, "TO": str(T), "TI": str(T), "WMIN": str(WMIN),
               "LMIN": str(LMIN), "workers": WORKERS, "chunk": CHUNK,
               "outer_nodes": N, "chunks": nch, "elapsed_s": el,
               "d0_err": mp.nstr(abs(mp.re(S[0])-S1), 3)},
              open("modelA_h768_meta.json", "w"), indent=1)
    log("jobA768 DONE")

if __name__ == "__main__":
    main()
'''
open("jobA_h768.py","w").write(job)
import py_compile; py_compile.compile("jobA_h768.py", doraise=True); print("jobA_h768.py written, compiles OK")

step06 = r'''"""step_06.py -- launch jobA_h768.py DETACHED (survives this session) at nice 10.

Command (exact, also in result.md):
  nohup nice -n 10 $VENV/python3 $TOOLS/detach.py $VENV/python3 \
      $CWD/jobA_h768.py > $CWD/jobA_h768_out.txt 2>&1 &
detach.py does os.setsid()+execvp, so the job gets its own session; the nice
level propagates across exec and is inherited by the 6 spawn workers.
"""
import subprocess, sys, os, time, json
CWD = "<project>/research/08_theory/K4h_tetrakis_parallel"
VENV = "<project>/.venv/bin/python3"
DETACH = "<project>/tools/detach.py"
JOB = os.path.join(CWD, "jobA_h768.py")

def main():
    assert os.path.exists(VENV), VENV
    assert os.path.exists(DETACH), DETACH
    assert os.path.exists(JOB), JOB
    out = open(os.path.join(CWD, "jobA_h768_out.txt"), "w")
    cmd = ["nohup", "nice", "-n", "10", VENV, DETACH, VENV, JOB]
    p = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, cwd=CWD)
    print("launcher PID (detach.py, becomes the job after exec):", p.pid)
    print("command:", " ".join(cmd), "&")
    time.sleep(5)
    pidf = os.path.join(CWD, "jobA_h768.pid")
    if os.path.exists(pidf):
        print("job pid file:", open(pidf).read())
    json.dump({"launcher_pid": p.pid, "cmd": " ".join(cmd) + " &"},
              open(os.path.join(CWD, "launch_info.json"), "w"))

if __name__ == "__main__":
    main()
'''
open("step_06.py","w").write(step06)
import py_compile; py_compile.compile("step_06.py", doraise=True); print("step_06.py written, compiles OK")
import time; time.sleep(120)
print(open("val96_out.txt").read())
