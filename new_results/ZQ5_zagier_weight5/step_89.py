import sys, subprocess
lib = open('w5lib.py').read()
lib = lib.replace("""def int_left_kernel(A):""", """def int_left_kernel_hnf(A):
    \"\"\"Saturated Z-basis of {c : c A = 0} via HNF with transformation (fast, flint).\"\"\"
    import flint as _fl
    M = _fl.fmpz_mat([list(map(int, r)) for r in A])
    H, U = M.hnf(transform=True)
    nr = M.nrows(); nc = M.ncols()
    rank = 0
    for i in range(nr):
        if any(H[i, j] for j in range(nc)): rank += 1
    K = []
    for i in range(rank, nr):
        K.append([int(U[i, j]) for j in range(nc if False else U.ncols())])
    # U is (nr x nr); rows rank..nr-1 satisfy u*M == 0
    return K

def int_left_kernel(A):""")
open('w5lib.py','w').write(lib)
tst = r'''
import sys; sys.path.insert(0,'.')
import sympy as sp
import w5lib as W
# test vs old method on a planted matrix
A = [[1,0],[2,1],[3,2],[1,1]]
K1 = W.int_left_kernel(A)
K2 = W.int_left_kernel_hnf(A)
print("old:", K1)
print("new:", K2)
assert sp.Matrix(K1).rank() == sp.Matrix(K2).rank()
assert sp.Matrix(K2)*sp.Matrix(A) == sp.zeros(len(K2), 2)
print("snf new:", W.snf_invariants(K2))
# D=8 kernel consistency
import pipe, json
C = W.letters_info(8)
J, cands = pipe.load_cands(8)
d = 2
A8 = [W.symbol_li5_dense(c['vx'], c['v1m'], d) for c in cands]
P,_ = W.sym5_basis(d)
Wann_old = W.int_left_kernel([list(col) for col in zip(*P)])
Wann_new = W.int_left_kernel_hnf([list(col) for col in zip(*P)])
print("ann dims old/new:", len(Wann_old), len(Wann_new))
Am = sp.Matrix(A8)
Apr_old = Am*sp.Matrix(Wann_old).T
Apr_new = Am*sp.Matrix(Wann_new).T
K_old = W.int_left_kernel([list(r) for r in Apr_old.tolist()])
K_new = W.int_left_kernel_hnf([list(r) for r in Apr_new.tolist()])
print("D=8 kernel dims old/new:", len(K_old), len(K_new))
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',tst],capture_output=True,text=True,timeout=900)
print(r.stdout); print("STDERR:", r.stderr[-1000:])
with open('progress.log','a') as f: f.write("step_54_hnfkernel: HNF kernel implemented+tested\n")
