
"""K4 engine: surface area of K_p = {sum_j |<n_j,x>|^p <= 1} in R^3 near p=1."""
import mpmath as mp
from mpmath import mpf

KMAX = 6

def ser_mul(a, b, K):
    out = [mpf(0)]*(K+1)
    for i in range(K+1):
        s = mpf(0)
        for j in range(i+1):
            s += a[j]*b[i-j]
        out[i] = s
    return out

def ser_inv(a, K):
    out = [mpf(0)]*(K+1)
    out[0] = 1/a[0]
    for n in range(1, K+1):
        s = mpf(0)
        for j in range(1, n+1):
            s += a[j]*out[n-j]
        out[n] = -s/a[0]
    return out

def ser_log(a, K):
    l0 = mp.log(a[0])
    b = [a[i]/a[0] for i in range(K+1)]
    u = b[:]; u[0] = mpf(0)
    out = [mpf(0)]*(K+1)
    pwr = [mpf(0)]*(K+1); pwr[0] = mpf(1)
    for m in range(1, K+1):
        pwr = ser_mul(pwr, u, K)
        cm = mpf(1 if m%2 else -1)/m   # (-1)^(m+1)/m
        for i in range(m, K+1):
            out[i] += cm*pwr[i]
    out[0] = l0
    return out

def ser_exp(a, K):
    out = [mpf(0)]*(K+1)
    out[0] = mp.exp(a[0])
    for n in range(1, K+1):
        s = mpf(0)
        for j in range(1, n+1):
            s += j*a[j]*out[n-j]
        out[n] = s/n
    return out

def ser_sqrt(a, K):
    s0 = mp.sqrt(a[0])
    b = [a[i]/a[0] for i in range(K+1)]
    u = b[:]; u[0] = mpf(0)
    out = [mpf(0)]*(K+1)
    pwr = [mpf(0)]*(K+1); pwr[0] = mpf(1)
    bc = [mpf(1)]
    for m in range(1, K+1):
        bc.append(bc[-1]*(mpf(1)/2 - (m-1))/m)
    for m in range(1, K+1):
        pwr = ser_mul(pwr, u, K)
        for i in range(m, K+1):
            out[i] += bc[m]*pwr[i]
    out[0] = mpf(1)
    return [s0*out[i] for i in range(K+1)]

class Model:
    def __init__(self, normals, K=KMAX):
        self.n = [mp.matrix(v) for v in normals]
        self.m = len(normals)
        self.K = K

    def c_vals(self, u):
        return [(self.n[j].T*u)[0] for j in range(self.m)]

    def F_series(self, u, sigmas):
        K = self.K
        c = self.c_vals(u)
        a = [sigmas[j]*c[j] for j in range(self.m)]
        L = [mp.log(a[j]) for j in range(self.m)]
        Phi = [mpf(0)]*(K+1)
        fk = [1/mp.factorial(k) for k in range(K+1)]
        for k in range(K+1):
            s = mpf(0)
            for j in range(self.m):
                s += a[j] * L[j]**k * fk[k]
            Phi[k] = s
        EE = []
        for k in range(K+1):
            v = mp.matrix(3,1)
            for j in range(self.m):
                coef = sigmas[j]*L[j]**k*fk[k]
                v += coef*(self.n[j] - c[j]*u)
            EE.append(v)
        Wdot = [mpf(0)]*(K+1)
        for k in range(K+1):
            s = mpf(0)
            for i in range(k+1):
                s += (EE[i].T*EE[k-i])[0]
            Wdot[k] = s
        Phi2inv = ser_inv(ser_mul(Phi, Phi, K), K)
        W = ser_mul(Wdot, Phi2inv, K)
        logPhi = ser_log(Phi, K)
        inv1pd = [mpf(-2)*(-1)**k for k in range(K+1)]
        expo = ser_mul(inv1pd, logPhi, K)
        factor = ser_exp(expo, K)
        oneW = W[:]; oneW[0] += 1
        sq = ser_sqrt(oneW, K)
        return ser_mul(factor, sq, K)

    def F_direct(self, u, sigmas, d):
        c = self.c_vals(u)
        a = [sigmas[j]*c[j] for j in range(self.m)]
        p = 1+d
        Phi = sum(a[j]**p for j in range(self.m))
        B = mp.matrix(3,1)
        for j in range(self.m):
            B += (a[j]**d)*sigmas[j]*(self.n[j] - c[j]*u)
        V = B/Phi
        W = (V.T*V)[0]
        return Phi**(-2/p)*mp.sqrt(1+W)

def tetra_model(K=KMAX):
    s3 = 1/mp.sqrt(3)
    return Model([[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]], K)

def octa_model(K=KMAX):
    return Model([[1,0,0],[0,1,0],[0,0,1]], K)

def th_A(ph):
    return mp.atan(1/mp.sin(ph))
def th_B(ph):
    return mp.atan(1/(mp.cos(ph)+mp.sin(ph)))

SIG_A = [+1,+1,+1,+1]
SIG_B = [+1,-1,+1,+1]
SIG_OCT = [+1,+1,+1]


# ---------- fast versions (plain tuples, no mp.matrix) ----------
class FastModel:
    def __init__(self, normals, K=4):
        self.n = [tuple(mpf(x) for x in v) for v in normals]
        self.m = len(normals)
        self.K = K
        G = [[mpf(0)]*self.m for _ in range(self.m)]
        for i in range(self.m):
            for j in range(self.m):
                G[i][j] = sum(self.n[i][t]*self.n[j][t] for t in range(3))
        self.G = G
        self.fk = [1/mp.factorial(k) for k in range(K+1)]
        # binomial sqrt coeffs
        bc = [mpf(1)]
        for m in range(1, K+1):
            bc.append(bc[-1]*(mpf(1)/2 - (m-1))/m)
        self.bc = bc

    def F_series(self, ux, uy, uz, sigmas):
        K = self.K; m = self.m; fk = self.fk
        n = self.n
        c = [n[j][0]*ux + n[j][1]*uy + n[j][2]*uz for j in range(m)]
        a = [sigmas[j]*c[j] for j in range(m)]
        L = [mp.log(a[j]) for j in range(m)]
        # Phi series
        Phi = [mpf(0)]*(K+1)
        for k in range(K+1):
            s = mpf(0)
            for j in range(m):
                s += a[j]*L[j]**k*fk[k]
            Phi[k] = s
        # EE[k] vector (3 comps)
        EEx = [mpf(0)]*(K+1); EEy = [mpf(0)]*(K+1); EEz = [mpf(0)]*(K+1)
        for k in range(K+1):
            sx=sy=sz=mpf(0)
            for j in range(m):
                coef = sigmas[j]*L[j]**k*fk[k]
                bx = n[j][0]-c[j]*ux; by = n[j][1]-c[j]*uy; bz = n[j][2]-c[j]*uz
                sx += coef*bx; sy += coef*by; sz += coef*bz
            EEx[k]=sx; EEy[k]=sy; EEz[k]=sz
        Wdot = [mpf(0)]*(K+1)
        for k in range(K+1):
            s = mpf(0)
            for i in range(k+1):
                s += EEx[i]*EEx[k-i]+EEy[i]*EEy[k-i]+EEz[i]*EEz[k-i]
            Wdot[k] = s
        P2 = ser_mul(Phi, Phi, K)
        P2i = ser_inv(P2, K)
        W = ser_mul(Wdot, P2i, K)
        logPhi = ser_log(Phi, K)
        inv1pd = [mpf(-2)*(-1)**k for k in range(K+1)]
        expo = ser_mul(inv1pd, logPhi, K)
        factor = ser_exp(expo, K)
        oneW = W[:]; oneW[0] += 1
        sq = ser_sqrt(oneW, K)
        return ser_mul(factor, sq, K)

    def F_direct(self, ux, uy, uz, sigmas, d):
        n = self.n; m = self.m
        c = [n[j][0]*ux + n[j][1]*uy + n[j][2]*uz for j in range(m)]
        a = [sigmas[j]*c[j] for j in range(m)]
        p = 1+d
        Phi = sum(a[j]**p for j in range(m))
        Vx=Vy=Vz=mpf(0)
        for j in range(m):
            coef = (a[j]**d)*sigmas[j]
            Vx += coef*(n[j][0]-c[j]*ux); Vy += coef*(n[j][1]-c[j]*uy); Vz += coef*(n[j][2]-c[j]*uz)
        W = (Vx*Vx+Vy*Vy+Vz*Vz)/(Phi*Phi)
        return Phi**(-2/p)*mp.sqrt(1+W)

def fast_tetra(K=4):
    s3 = 1/mp.sqrt(3)
    return FastModel([[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]], K)
def fast_octa(K=4):
    return FastModel([[1,0,0],[0,1,0],[0,0,1]], K)
