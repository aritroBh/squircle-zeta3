"""lllkernel.py -- approximate integer kernel of a real matrix via LLL (python-flint).
Given A (M x N) mpf rows, K scale: returns short integer vectors q with A q ~ 0."""
import mpmath as mp
from flint import fmpz_mat, fmpz

def int_kernel(A, Kexp=45):
    """A: list of rows (each a list of mpf), M x N. Returns list of integer vectors q (len N)."""
    M = len(A); N = len(A[0])
    K = mp.power(10, Kexp)
    # B = [I_N | round(K * A^T)] : N rows, N+M cols
    B = fmpz_mat(N, N + M)
    for i in range(N):
        B[i, i] = fmpz(1)
    for j in range(M):
        for i in range(N):
            x = mp.nint(K*A[j][i])
            B[i, N + j] = fmpz(int(x))
    BL = B.lll()
    out = []
    for i in range(N):
        q = [int(BL[i, j]) for j in range(N)]
        tail = [int(BL[i, N + j]) for j in range(M)]
        out.append((q, max(abs(t) for t in tail) if tail else 0))
    return out
