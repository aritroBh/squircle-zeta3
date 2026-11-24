step03 = r'''"""step_03.py -- cheap-params three-way equivalence test for k4h_driver.

(a) str(mpf) round-trip at dps 80 and 270 (checkpoint format safety).
(b) build_outer_nodes determinism.
(c) h=1/24, dps 50, WMIN 1e-30, LMIN 1e-20, T=4.0:
    original k4c_driver.chamber_S  vs  sequential node loop (inner_sum)  vs
    parallel Pool(3) with chunk=32 + combine.  Require max diff << 1e-45.
"""
import sys, os, time
sys.dont_write_bytecode = True
os.chdir("<project>/research/08_theory/K4h_tetrakis_parallel")
K4E = "<project>/research/08_theory/K4e_tetrakis_d3"
sys.path.insert(0, K4E)
import mpmath as mp
from mpmath import mpf

# (a) str round-trip
for dps in (80, 270):
    mp.mp.dps = dps
    ok = True
    x = mpf(1)
    for i in range(200):
        x = mp.sin(x + mp.pi/7) + mp.log(x+2) + mp.sqrt(x+1)
        if mpf(str(x)) != x: ok = False
    print("str(mpf) round-trip dps=%d: %s" % (dps, "EXACT" if ok else "FAIL"))

mp.mp.dps = 50
import k4c_driver
from k4c_driver import chamber_S
import k4h_driver as KD
import multiprocessing as M

s2 = mp.sqrt(2)
NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]

# (b) determinism of node list
n1 = KD.build_outer_nodes(NA, mpf(1)/24, mpf('4.0'), mpf('1e-30'))
n2 = KD.build_outer_nodes(NA, mpf(1)/24, mpf('4.0'), mpf('1e-30'))
det = len(n1) == len(n2) and all(str(a[0])==str(b[0]) and str(a[1])==str(b[1]) for a,b in zip(n1,n2))
print("node list determinism: %s (N=%d)" % ("OK" if det else "FAIL", len(n1)))

# (c) three-way equivalence
WMIN = mpf('1e-30'); LMIN = mpf('1e-20')
k4c_driver.WMIN = WMIN; k4c_driver.LMIN = LMIN
t0 = time.time()
S_ref, st_ref = chamber_S(NA, 4, mpf(1)/24, mpf(1)/24, mpf('4.0'), mpf('4.0'), HI2=mpf(1)/24)
t_ref = time.time()-t0
print("original chamber_S done: %.1fs stats=%s" % (t_ref, st_ref))

# sequential node loop through inner_sum
t0 = time.time()
Mm = KD.FastModel(NA, 4)
acc = [mpf(0)]*5
st_seq = {'zero_len':0,'nbp':0,'nsing':0,'ninner':0}
for (ph,w) in n1:
    a, st = KD.inner_sum(Mm, ph, w, 4, mpf(1)/24, mpf('4.0'), mpf(1)/24, WMIN, LMIN)
    for k in range(5): acc[k] += a[k]
    for key in st_seq: st_seq[key] += st[key]
S_seq = [48*v for v in acc]
print("sequential node loop done: %.1fs stats=%s" % (time.time()-t0, st_seq))

# parallel
cfg = dict(normals=[[str(x) for x in row] for row in NA], K=4,
           HO=str(mpf(1)/24), HI=str(mpf(1)/24), HI2=str(mpf(1)/24),
           TO='4.0', TI='4.0', WMIN=str(WMIN), LMIN=str(LMIN), dps=50, chunk=32)
nch = (len(n1)+31)//32
t0 = time.time()
ctx = M.get_context("spawn")
with ctx.Pool(3, initializer=KD.init_worker, initargs=(cfg,)) as pool:
    res = list(pool.imap_unordered(KD.do_chunk, range(nch)))
acc3, st3 = KD.combine(res, 4)
S_par = [48*v for v in acc3]
print("parallel (3 workers) done: %.1fs stats=%s" % (time.time()-t0, st3))

print()
for k in range(5):
    d1 = abs(mp.re(S_seq[k]) - mp.re(S_ref[k]))
    d2 = abs(mp.re(S_par[k]) - mp.re(S_ref[k]))
    print("k=%d  |seq-orig|=%s   |par-orig|=%s" % (k, mp.nstr(d1,3), mp.nstr(d2,3)))
print("stats match (3 shared keys):",
      all(st_seq[key]==st_ref[key] for key in ('zero_len','nbp','nsing')),
      all(st3[key]==st_ref[key] for key in ('zero_len','nbp','nsing')))
'''
open("step_03.py","w").write(step03)
import subprocess, sys
r = subprocess.run([sys.executable, "step_03.py"], capture_output=True, text=True)
print(r.stdout); print("ERR:", r.stderr[-3000:] if r.returncode else "(clean)")
