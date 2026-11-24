import flint, mpmath, json, time
from mpmath import mp, mpf, pi, nstr
print("flint version ok")
# test fmpz_mat lll
M = flint.fmpz_mat([[1,2],[3,4]])
print(hasattr(M, 'lll'))
import inspect
try:
    print(M.lll())
except Exception as e:
    print("lll err:", e)
