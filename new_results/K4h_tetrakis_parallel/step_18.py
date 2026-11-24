step05 = r'''"""step_05.py -- K4h task 3: production-cost pilot at full production params.

Params: h=1/768, dps 270, WMIN 1e-215, LMIN 1e-205, TO=TI=6.0 (step_02:
t*=5.75 for these WMIN/h; w(6.0)=4e-276 << WMIN, so the weights reach WMIN
before T cuts -- the requested 5.8-6.0 window confirmed).

Measures per-outer-node cost (time + inner-eval count) on a strided sample of
60 nodes through a 6-worker spawn Pool, then extrapolates to all N nodes.
"""
import sys, os, time, json
sys.dont_write_bytecode = True
os.chdir("<project>/research/08_theory/K4h_tetrakis_parallel")
K4E = "<project>/research/08_theory/K4e_tetrakis_d3"
if K4E not in sys.path: sys.path.insert(0, K4E)
import mpmath as mp
from mpmath import mpf
import multiprocessing as M

def main():
    mp.mp.dps = 270
    import k4h_driver as KD
    s2 = mp.sqrt(2)
    NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
    cfg = dict(normals=[[str(x) for x in row] for row in NA], K=4,
               HO=str(mpf(1)/768), HI=str(mpf(1)/768), HI2=str(mpf(1)/768),
               TO='6.0', TI='6.0', WMIN='1e-215', LMIN='1e-205', dps=270, chunk=24)
    t0 = time.time()
    nodes = KD.build_outer_nodes(NA, mpf(1)/768, mpf('6.0'), mpf('1e-215'))
    N = len(nodes)
    print("outer nodes N=%d built in %.1fs" % (N, time.time()-t0), flush=True)
    NS = 60
    sample = [round(i*(N-1)/(NS-1)) for i in range(NS)]
    ctx = M.get_context("spawn")
    t0 = time.time()
    with ctx.Pool(6, initializer=KD.init_worker, initargs=(cfg,)) as pool:
        res = list(pool.imap_unordered(KD.do_nodes, [[i] for i in sample]))
    wall = time.time()-t0
    els = sorted(r['el'] for r in res)
    nin = sorted(r['stats']['ninner'] for r in res)
    nzero = sum(1 for r in res if r['stats']['zero_len'] and not r['stats']['ninner'])
    tot_el = sum(els); tot_nin = sum(nin)
    print("sample: %d nodes, wall=%.0fs with 6 workers" % (NS, wall))
    print("per-node time s: min=%.1f med=%.1f mean=%.1f max=%.1f" %
          (els[0], els[NS//2], tot_el/NS, els[-1]))
    print("inner evals per node: min=%d med=%d mean=%.0f max=%d ; zero_len-only nodes: %d" %
          (nin[0], nin[NS//2], tot_nin/NS, nin[-1], nzero))
    mean_node = tot_el/NS
    single_total = mean_node*N
    print("ESTIMATE single-process total: %.0f s = %.1f h" % (single_total, single_total/3600))
    print("ESTIMATE 6-worker wall (mean-based): %.0f s = %.1f h" % (single_total/6, single_total/3600/6))
    wall_est = wall*N/NS
    print("ESTIMATE 6-worker wall (throughput-based): %.0f s = %.1f h" % (wall_est, wall_est/3600))
    # chunks of 24: checkpoint cadence
    nch = (N+23)//24
    per_chunk = mean_node*24
    print("chunks=%d (size 24), per-chunk ~%.0f s (~%.1f min); checkpoint cadence ~%.1f min"
          % (nch, per_chunk, per_chunk/60, per_chunk/60/6))
    json.dump({"N": N, "sample": NS, "wall": wall, "mean_node": mean_node,
               "median_node": els[NS//2], "max_node": els[-1], "mean_ninner": tot_nin/NS,
               "est_single_s": single_total, "est_6w_mean_s": single_total/6,
               "est_6w_throughput_s": wall_est, "nch24": nch},
              open("pilot_h768.json", "w"), indent=1)
    print("saved pilot_h768.json")

if __name__ == "__main__":
    main()
'''
open("step_05.py","w").write(step05)
import py_compile; py_compile.compile("step_05.py", doraise=True); print("step_05 written, compiles OK")
import time
time.sleep(60)
print(open("val96_out.txt").read())
