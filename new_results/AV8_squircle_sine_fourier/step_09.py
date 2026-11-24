from mpmath import mp, mpf, pi, matrix, lu_solve, factorial
import json
mp.dps = 45
data = json.load(open('scratchpad/bn_near2.json'))
h = mpf('0.02')
offsets = list(range(-4,5))
ns = [1,3,5,7,9,11,13]

def derivs(n):
    # fit p(e) = sum c_k e^k through 9 points; c_k = b^{(k)}/k!
    A = matrix(9,9); b = matrix(9,1)
    for i,j in enumerate(offsets):
        for k in range(9):
            A[i,k] = (j*h)**k
        b[i] = mpf(data[f"{n},{j}"])
    c = lu_solve(A,b)
    return [c[k]*factorial(k) for k in range(9)]

res = {}
print(f"{'n':>3} {'b_n^(1)':>28} {'b_n^(2)':>28} {'b_n^(3)':>28}")
for n in ns:
    d = derivs(n)
    res[n] = d
    print(f"{n:>3} {mp.nstr(d[1],26):>28} {mp.nstr(d[2],26):>28} {mp.nstr(d[3],26):>28}")
print()
print("b_3'''/6 =", mp.nstr(res[3][3]/6, 38))
print("target   = -0.09771983272784923994573")
# rough error estimate from 8th/9th deriv magnitude
for n in [1,3]:
    d = res[n]
    print(n, [mp.nstr(d[k],5) for k in range(4,9)])
