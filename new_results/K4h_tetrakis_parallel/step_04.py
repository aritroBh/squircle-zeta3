"""step_04.py -- K4h task 2: validate parallel driver vs single-process chamber_S.

Params: h=1/96, dps 80, WMIN 1e-45, LMIN 1e-35 (the jobA_h384 convention:
WMIN/LMIN = 1e-(dps-35)/1e-(dps-45)), TO=TI=4.4 (step_02: t*=4.19 for these
WMIN/h, so weights reach WMIN well before T cuts: w(4.4)=2.3e-56).

Runs: (1) parallel run_job with 6 workers, chunk=32, checkpoints in ckpt_val96/;
(2) resume test: delete last chunk, rerun run_job -> must reproduce identical S;
(3) original single-process k4c_driver.chamber_S at identical params.
PASS iff max_k |S_par - S_ref| <= 1e-75.
"""
import sys, os, time, json, glob
sys.dont_write_bytecode = True
os.chdir("<project>/research/08_theory/K4h_tetrakis_parallel")
K4E = "<project>/research/08_theory/K4e_tetrakis_d3"
if K4E not in sys.path: sys.path.insert(0, K4E)
import mpmath as mp
from mpmath import mpf

LOGF = "val96.log"
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    open(LOGF, "a").write(s + "\n")

def main():
    mp.mp.dps = 80
    import k4c_driver
    from k4c_driver import chamber_S
    import k4h_driver as KD

    s2 = mp.sqrt(2)
    NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
    WMIN = mpf('1e-45'); LMIN = mpf('1e-35'); T = mpf('4.4'); H = mpf(1)/96
    cfg = dict(normals=[[str(x) for x in row] for row in NA], K=4,
               HO=str(H), HI=str(H), HI2=str(H), TO=str(T), TI=str(T),
               WMIN=str(WMIN), LMIN=str(LMIN), dps=80, chunk=32)
    open(LOGF, "w").write("val96 start dps=80 h=1/96 T=4.4 WMIN=1e-45 LMIN=1e-35 chunk=32 workers=6\n")
    log("[val96] cfg:", {k: cfg[k] for k in ('HO','TO','WMIN','LMIN','dps','chunk')})

    t0 = time.time()
    acc, stats, N, nch = KD.run_job(cfg, "ckpt_val96", 6, log, "val96")
    S_par = [48*v for v in acc]
    log("[val96] parallel done: %.0fs N=%d nch=%d stats=%s" % (time.time()-t0, N, nch, stats))
    for k in range(5):
        log("[val96] [d^%d]S_par = %s" % (k, mp.nstr(mp.re(S_par[k]), 70)))

    # ---- resume test: delete last chunk file, rerun, must give identical S ----
    files = sorted(glob.glob("ckpt_val96/chunk_*.json"))
    victim = files[-1]
    os.remove(victim)
    log("[val96-resume] deleted %s; rerunning run_job" % victim)
    t0 = time.time()
    acc2, stats2, _, _ = KD.run_job(cfg, "ckpt_val96", 6, log, "val96-resume")
    S_par2 = [48*v for v in acc2]
    dres = max(abs(mp.re(S_par2[k]) - mp.re(S_par[k])) for k in range(5))
    log("[val96-resume] rerun done: %.0fs; max |S_resume - S_par| = %s"
        % (time.time()-t0, mp.nstr(dres, 3)))

    # ---- single-process reference ----
    k4c_driver.WMIN = WMIN; k4c_driver.LMIN = LMIN
    t0 = time.time()
    S_ref, st_ref = chamber_S(NA, 4, H, H, T, T, HI2=H)
    log("[val96] sequential chamber_S done: %.0fs stats=%s" % (time.time()-t0, st_ref))
    for k in range(5):
        log("[val96] [d^%d]S_ref = %s" % (k, mp.nstr(mp.re(S_ref[k]), 70)))

    log("[val96] stats match:", all(stats[k]==st_ref[k] for k in ('zero_len','nbp','nsing')))
    worst = mpf(0)
    for k in range(5):
        d = abs(mp.re(S_par[k]) - mp.re(S_ref[k]))
        worst = max(worst, d)
        log("[val96] k=%d  |S_par - S_ref| = %s" % (k, mp.nstr(d, 3)))
    verdict = "PASS" if worst <= mpf('1e-75') else "FAIL"
    log("[val96] max diff = %s  -> %s (bound 1e-75)" % (mp.nstr(worst, 3), verdict))
    json.dump({"Spar": [mp.nstr(mp.re(v), 80) for v in S_par],
               "Sref": [mp.nstr(mp.re(v), 80) for v in S_ref],
               "maxdiff": mp.nstr(worst, 3), "verdict": verdict,
               "stats_par": stats, "stats_ref": st_ref},
              open("val96_result.json", "w"), indent=1)
    log("[val96] DONE")

if __name__ == "__main__":
    main()
