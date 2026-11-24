
import mpmath as mp
mp.mp.dps = 45
import json, sys

def S_p(phi, m, p):
    return mp.fsum(abs(mp.cos(phi - i*mp.pi/m))**p for i in range(1, m+1))
def rho_p(phi, m, p):
    return S_p(phi, m, p)**(-mp.mpf(1)/p)
def h_p(theta, m, p):
    f = lambda phi: rho_p(phi, m, p)*mp.cos(theta - phi)
    N = 400
    best = mp.mpf(-1); bph = 0
    # search only near theta
    for k in range(N):
        ph = (theta - mp.pi/2) + mp.pi*k/(N-1)
        v = f(ph)
        if v > best: best, bph = v, ph
    a, b = bph - mp.pi/N, bph + mp.pi/N
    gr = (mp.sqrt(5)-1)/2
    c, d = b - gr*(b-a), a + gr*(b-a)
    fc, fd = f(c), f(d)
    for _ in range(220):
        if fc < fd:
            a, c, fc = c, d, fd
            d = a + gr*(b-a); fd = f(d)
        else:
            b, d, fd = d, c, fc
            c = b - gr*(b-a); fc = f(c)
    return f((a+b)/2)
def tent(u, m):
    sm = mp.sin(mp.pi/(2*m))
    return sm*mp.cos(u) + sm*mp.tan(mp.pi/(2*m))*abs(mp.sin(u))

def Dval(m, p, s, tau):
    u = tau*s
    return mp.log(h_p(mp.pi/2 - mp.pi/(2*m) + u, m, p)/tent(u, m))/tau

def extract_profiles(m, s_list, tau0):
    # returns D1, D2, D3 dicts via Richardson in tau at each s
    out = []
    for s in s_list:
        vals = []
        for j in range(4):
            tau = tau0/2**j
            p = 1/(1-tau)
            vals.append(Dval(m, p, s, tau))
        # Richardson: D = D1 + tau D2 + tau^2 D3 + tau^3 D4
        D1 = mp.richardson if False else None
        # Neville at tau=0 with points (tau_j, vals_j)
        def nev(pts):
            n = len(pts); Q = [[mp.mpf(0)]*n for _ in range(n)]
            for i in range(n): Q[i][0] = pts[i][1]
            for jj in range(1,n):
                for i in range(n-jj):
                    xi, xj = pts[i][0], pts[i+jj][0]
                    Q[i][jj] = (xj*Q[i][jj-1] - xi*Q[i+1][jj-1])/(xj-xi)
            return Q[0][n-1]
        taus = [tau0/2**j for j in range(4)]
        D1 = nev(list(zip(taus, vals)))
        v2 = [(t, (v - D1)/t) for t, v in zip(taus, vals)]
        D2 = nev(v2)
        v3 = [(t, (v - D1 - D2*t)/t**2) for t, v in zip(taus, vals)]
        D3 = nev(v3)
        out.append((s, D1, D2, D3))
    return out

if __name__ == "__main__":
    mm = int(sys.argv[1])
    s_list = [mp.mpf(k)/4 for k in range(-80, 81)]  # s in [-20, 20] step 0.25
    prof = extract_profiles(mm, s_list, mp.mpf('0.006'))
    json.dump([[str(s), str(a), str(b), str(c)] for s,a,b,c in prof], open(f"ckpt_prof_m{mm}.json","w"))
    for s,a,b,c in prof[::20]:
        print(f"m={mm} s={mp.nstr(s,6)}: D1={mp.nstr(a,16)} D2={mp.nstr(b,16)} D3={mp.nstr(c,16)}")
