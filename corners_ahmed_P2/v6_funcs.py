import mpmath as mp

def dLam_fast(m, c, P, H):
    """Lambda'_m via Lambda_nu = (P_{nu-1} - P_{nu+1})/(2nu+1) and the degree derivative at integers
       d/dnu P_nu |_{nu=n} = P_n(c)[log((1+c)/2) + 2(H_{2n} - H_n)] + 2 sum_{k<n} (-1)^{n+k}(2k+1)P_k(c)/((n-k)(n+k+1))
       (Szmytkowski 2006; re-verified against numerical differentiation in V6e).  P, H: precomputed tables."""
    def dP(n):
        if n < 0:            # P_mu = P_{-mu-1}  =>  dP at mu = n equals -dP at -n-1
            return -dP(-n - 1)
        s = P[n]*(mp.log((1 + c)/2) + 2*(H[2*n] - H[n]))
        s += 2*mp.fsum((-1)**(n + k)*(2*k + 1)*P[k]/((n - k)*(n + k + 1)) for k in range(n))
        return s
    Lm = (P[m - 1] - P[m + 1])/(2*m + 1) if m >= 1 else 1 - c
    return (dP(m - 1) - dP(m + 1))/(2*m + 1) - 2*Lm/(2*m + 1)


def tables(c, N):
    P = [mp.mpf(1), c]
    for k in range(1, N + 1):
        P.append(((2*k + 1)*c*P[k] - k*P[k - 1])/(k + 1))
    H = [mp.mpf(0)]
    for k in range(1, 2*N + 4):
        H.append(H[-1] + mp.mpf(1)/k)
    return P, H


def Kalpha(th):
    k = mp.sin(th/2)
    return -(1 + k)*mp.log(1 + k) - (1 - k)*mp.log(1 - k)


