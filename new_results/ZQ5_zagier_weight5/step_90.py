import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=17
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
d = len(C['letters'])
N = len(cands)
t0=time.time()
A = [W.symbol_li5_dense(c['vx'], c['v1m'], d) for c in cands]
P, mskeys = W.sym5_basis(d)
Wann = W.int_left_kernel_hnf([list(col) for col in zip(*P)])
import sympy as sp
Am = sp.Matrix(A); Wm = sp.Matrix(Wann)
Apr = Am*Wm.T
t1=time.time()
print(f"symbol setup {t1-t0:.0f}s: tensor {d**5}, Sym {len(P)}, quotient {len(Wann)}")
K = W.int_left_kernel_hnf([list(r) for r in Apr.tolist()])
Km = sp.Matrix(K)
assert Km*Apr == sp.zeros(len(K), Apr.shape[1])
print(f"kernel dim {len(K)}, rank(Apr)={Apr.rank()}, snf={W.snf_invariants(K)[:5]}... ({time.time()-t1:.0f}s)")
json.dump(dict(K=K), open(pipe.cand_dir(D)+"/kernel.json","w"))
t0=time.time()
vals = pipe.phase_numerics(D, C, cands, P=1150)
print(f"numerics {time.time()-t0:.0f}s")
mp.mp.dps=60
print("L(chi17,5) =", mp.nstr(vals["L"], 30))
print("T = L/sqrt17 =", mp.nstr(vals["T"], 30))
'''
open('steps/step_54_d17_kn.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_54_d17_kn.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout); print("STDERR:", r.stderr[-1200:])
with open('progress.log','a') as f: f.write("step_54_d17_kn: D=17 kernel(HNF)+numerics\n")
