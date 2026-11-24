import time; t0=time.time()
import flint
from flint import fmpq_mat, fmpq
# API probe
M = fmpq_mat(2,3,[1,2,3,4,5,6])
print('rref:', M.rref())
try:
    print('nullspace:', M.nullspace())
except Exception as e:
    print('nullspace err:', e)
A = fmpq_mat(2,2,[2,1,1,1])
b = fmpq_mat(2,1,[1,1])
print('solve:', A.solve(b))
print('rank:', M.rank())
print('t=%.1f'%(time.time()-t0))
