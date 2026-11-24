import os, json, time
os.chdir("research/08_theory/K4b_cubo_finish")

# ---- step_04.py: patch k4_quad.py with weight-thresholded node skipping (fix part 1) ----
src = open("k4_quad.py").read()
if "wmin" not in src:
    src = src.replace(
"""def tanh_sinh_nodes(a, b, h, T):
    \"\"\"nodes (x,w) for int_a^b f dx; t = k h, |t|<=T.\"\"\"""",
"""def tanh_sinh_nodes(a, b, h, T, wmin=None):
    \"\"\"nodes (x,w) for int_a^b f dx; t = k h, |t|<=T.
    K4b fix: nodes with |w| <= wmin are skipped (their total mass is
    negligible and their x-coordinates can carry no valid digits).\"\"\"""")
    src = src.replace(
"""    out = []
    for (x,w,k) in nodes:
        if k == 0:
            out.append((x,w))
        else:
            if x > a and x < b:
                out.append((x,w))
            x2 = a+b-x
            if x2 > a and x2 < b:
                out.append((x2, w))
    return out""",
"""    out = []
    for (x,w,k) in nodes:
        if wmin is not None and abs(w) <= wmin:
            continue
        if k == 0:
            out.append((x,w))
        else:
            if x > a and x < b:
                out.append((x,w))
            x2 = a+b-x
            if x2 > a and x2 < b:
                out.append((x2, w))
    return out""")
    open("k4_quad.py","w").write(src)
print(open("k4_quad.py").read())

import mpmath as mp
from mpmath import mpf
mp.mp.dps = 90
import importlib, k4_engine, k4_quad
importlib.reload(k4_engine); importlib.reload(k4_quad)
from k4_engine import *
from k4_quad import *

# timing probe
MT = fast_tetra(K=4)
ph, th = mpf('0.9'), mpf('0.5')
u = (mp.sin(th)*mp.cos(ph), mp.sin(th)*mp.sin(ph), mp.cos(th))
t0=time.time()
for _ in range(50): MT.F_series(u[0],u[1],u[2], SIG_A)
t_ser = (time.time()-t0)/50
t0=time.time()
for _ in range(50): MT.F_direct(u[0],u[1],u[2], SIG_A, mpf('0.01'))
t_dir = (time.time()-t0)/50
n_out = len(tanh_sinh_nodes(mp.pi/4, mp.pi/2, mpf(1)/32, mpf('4.5')))
n_in  = len(tanh_sinh_nodes(mpf(0), mpf('0.6'), mpf(1)/32, mpf('4.5')))
print(f"dps=90  F_series: {t_ser*1e3:.2f} ms   F_direct: {t_dir*1e3:.3f} ms")
print(f"nodes outer={n_out} inner~{n_in}  -> est. full run {n_out*n_in*2*t_ser/60:.1f} min (series), FD point {n_out*n_in*2*t_dir:.1f} s")
open("step_04.py","w").write(open(__file__).read() if "__file__" in dir() else "# see transcript\n")
