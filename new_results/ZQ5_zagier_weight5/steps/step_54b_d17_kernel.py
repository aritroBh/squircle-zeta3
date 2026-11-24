
import sys; sys.path.insert(0,'.')
import json, time
import w5lib as W
import pipe
import flint
D=17
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
d = len(C['letters']); N = len(cands)
t0=time.time()
A = [W.symbol_li5_dense(c['vx'], c['v1m'], d) for c in cands]
P, mskeys = W.sym5_basis(d)
Wann = W.int_left_kernel_hnf([list(col) for col in zip(*P)])
print(f"ann {time.time()-t0:.0f}s dim {len(Wann)}")
t0=time.time()
Af = flint.fmpz_mat(A)          # N x 1024
Wf = flint.fmpz_mat([list(r) for r in zip(*Wann)])  # 1024 x len(Wann)
Aprf = Af*Wf                    # N x len(Wann)
Apr = [[Aprf[i,j] for j in range(Aprf.ncols())] for i in range(Aprf.nrows())]
print(f"matmul {time.time()-t0:.0f}s")
t0=time.time()
K = W.int_left_kernel_hnf(Apr)
print(f"kernel {time.time()-t0:.0f}s dim {len(K)}")
# verify with flint
Kf = flint.fmpz_mat(K) if K else None
if K:
    chk = Kf*Aprf
    print("kernel check max|entry|:", max(abs(chk[i,j]) for i in range(chk.nrows()) for j in range(chk.ncols())))
print("snf:", W.snf_invariants(K)[:10])
json.dump(dict(K=K), open(pipe.cand_dir(D)+"/kernel.json","w"))
print("saved")
