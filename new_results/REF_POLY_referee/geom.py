import mpmath as mp
from jet import vertices
def edges(U):
    V = vertices(U); n = len(V); out = []
    for i in range(n):
        A, B = V[i], V[(i+1) % n]
        out.append((A, B))
    return out
def kink_of(U, V):
    vals = [abs(u[0]*V[0]+u[1]*V[1]) for u in U]
    j = min(range(len(U)), key=lambda i: vals[i]); assert vals[j] < mp.mpf(10)**(-mp.mp.dps//2)
    return j
def edge_data(U, A, B):
    L = mp.sqrt((B[0]-A[0])**2+(B[1]-A[1])**2)
    tau = ((B[0]-A[0])/L, (B[1]-A[1])/L)
    L0 = -(A[0]*tau[0]+A[1]*tau[1]); L1 = L-L0
    h = abs(A[0]*tau[1]-A[1]*tau[0])
    j0 = kink_of(U, A); j1 = kink_of(U, B)
    xi0 = abs(U[j0][0]*B[0]+U[j0][1]*B[1]); xi1 = abs(U[j1][0]*A[0]+U[j1][1]*A[1])
    return L, h, L0, L1, xi0, xi1, j0, j1
def Z_e(U, A, B):
    L, h, L0, L1, xi0, xi1, _, _ = edge_data(U, A, B)
    return 6*xi0*xi1*h**2/L*(3*(xi0*L0+xi1*L1)/L-2)
def sigma(U, V):
    s = 0
    for u in U:
        a = abs(u[0]*V[0]+u[1]*V[1])
        if a > mp.mpf(10)**(-mp.mp.dps//2): s += a*mp.log(a)
    return s
def P1_general(U):
    """referee's own first-order formula: P'(1) = sum_e sigma(V0)<tau,V0> - sigma(V1)<tau,V1>"""
    tot = 0
    for A, B in edges(U):
        L = mp.sqrt((B[0]-A[0])**2+(B[1]-A[1])**2); tau = ((B[0]-A[0])/L, (B[1]-A[1])/L)
        tot += sigma(U, A)*(tau[0]*A[0]+tau[1]*A[1]) - sigma(U, B)*(tau[0]*B[0]+tau[1]*B[1])
    return tot
