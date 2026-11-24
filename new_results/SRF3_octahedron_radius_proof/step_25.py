
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
step04 = r'''
"""
SRF3 step_04: decisive test for the singularity content of the Borel transform.
c_k = (a_k/k!) R^k, R=|z0|.  If B(z) has a sqrt branch at z0 (edge) only:
   c_k ~ k^be [C1 cos(th_e k + ph1)].
If also a vertex singularity at (log2+-i pi)/2 (arg th_v), add k^bv C2 cos(th_v k+ph2).
We (i) apply the edge-annihilator (1 - 2cos(th_e)E + E^2) and examine what remains;
(ii) least-squares fit both models and compare residuals.
"""
import mpmath as mp, json, os
mp.mp.dps = 80
base = "<project>"
h1 = json.load(open(os.path.join(base,"research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json")))
tau = [mp.mpf(x) for x in h1['eta']]
K = 80
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0)
th_e = float(mp.arg(z0))          # 1.787953
th_v = float(mp.pi - mp.arg(z0))  # 1.353640
b = [tau[k]/mp.factorial(k) for k in range(K+1)]
c = [float(b[k]*R**k) for k in range(K+1)]

# (i) edge annihilator
an = [c[k+2] - 2*mp.cos(th_e)*c[k+1] + c[k] for k in range(K-1)]
print("edge-annihilated |A_k| (should decay like k^{be-2} if pure edge):")
for k in [20,30,40,50,60,70]:
    print(f"  k={k}: {abs(an[k]):.3e}   |c_k|={abs(c[k]):.3e}")
# fit residual to k^g cos(th_v k + ph): log |an| envelope
import math
xs=[k for k in range(30,K-2) if abs(an[k])>1e-12]
# envelope via local maxima
print("log|an| at sample k:", [(k, round(math.log10(abs(an[k])),2)) for k in range(40,75,5)])

# (ii) LSQ fits
import numpy as np
def design(ks, comps):
    # comps: list of (beta, theta). columns: k^beta cos(th k), k^beta sin(th k)
    M = np.zeros((len(ks), 2*len(comps)))
    for i,k in enumerate(ks):
        for j,(be,th) in enumerate(comps):
            M[i,2*j]   = k**be*math.cos(th*k)
            M[i,2*j+1] = k**be*math.sin(th*k)
    return M

def fit(ks, comps):
    M = design(ks, comps)
    y = np.array([c[k] for k in ks])
    x, res, rank, sv = np.linalg.lstsq(M, y, rcond=None)
    r = y - M@x
    return x, float(np.sqrt(np.mean(r*r))), float(np.max(np.abs(r)))

print("\n=== model A: edge only, beta free on a grid ===")
ks = list(range(30, K+1))
for be in [-2.0,-1.75,-1.5,-1.25,-1.0,-0.5]:
    x, rms, mx = fit(ks, [(be, th_e)])
    print(f"  beta={be}: rms={rms:.3e} max={mx:.3e}  C={math.hypot(x[0],x[1]):.4e}")
print("\n=== model B: edge+vertex, both beta=-1.5 ===")
x, rms, mx = fit(ks, [(-1.5, th_e),(-1.5, th_v)])
Ce = math.hypot(x[0],x[1]); Cv = math.hypot(x[2],x[3])
print(f"  rms={rms:.3e}  C_edge={Ce:.4e}  C_vertex={Cv:.4e}  ratio={Cv/Ce:.3e}")
print("\n=== model C: edge+vertex, beta grid ===")
for be in [-1.5,-1.0,-0.5]:
    for bv in [-1.5,-0.5,0.5]:
        x, rms, mx = fit(ks, [(be, th_e),(bv, th_v)])
        Ce = math.hypot(x[0],x[1]); Cv = math.hypot(x[2],x[3])
        print(f"  be={be} bv={bv}: rms={rms:.3e}  Ce={Ce:.3e} Cv={Cv:.3e}")
print("\n=== model A refined: joint (beta,theta) free near edge — via grid on theta ===")
best=None
for i in range(-20,21):
    th = th_e*(1+i*2e-6)
    x, rms, mx = fit(ks, [(-1.5, th)])
    if best is None or rms<best[0]: best=(rms,th,x)
print(f"  best theta={best[1]:.8f} (pred {th_e:.8f}) rms={best[0]:.3e}")
'''
open(wd+"/step_04.py","w").write(step04)
import subprocess
PY = "<project>/.venv/bin/python"
r = subprocess.run(["nice","-n","10",PY,"step_04.py"],capture_output=True,text=True,cwd=wd,timeout=900)
print(r.stdout); print(r.stderr[-1500:])
