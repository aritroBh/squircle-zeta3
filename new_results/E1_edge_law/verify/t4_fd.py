# T4: the task's suggested independent route: direct arc length / perimeter at p = 1 + k*h (no d-expansion at all),
# then a symmetric (2K+1)-point polynomial fit for P'''(1).  Rational hexagon b=0 (arc 1 = the r=2 edge, and the total).
import sys, json, time, mpmath as mp
from multiprocessing import Pool
sys.path.insert(0, 'research/10_frontier/E1_edge_law/verify')
DPS = int(sys.argv[1]); mp.mp.dps = DPS; H = mp.mpf(sys.argv[2]); K = int(sys.argv[3])
mp.mp.dps = DPS
vs = [(mp.mpf(0), mp.mpf(1)/4), (-mp.mpf(1)/4, mp.mpf(1)/4), (-mp.mpf(1)/4, mp.mpf(0))]
cs = [mp.mpf(1), mp.mpf(3), mp.mpf(1)]
def integrand(th, p):
    u = (mp.cos(th), mp.sin(th)); du = (-mp.sin(th), mp.cos(th))
    F = mp.mpf(0); dF = mp.mpf(0)
    for v, c in zip(vs, cs):
        g = v[0]*u[0] + v[1]*u[1]; gp = v[0]*du[0] + v[1]*du[1]
        s = abs(g); sg = 1 if g > 0 else -1
        F += c*s**p; dF += c*p*s**(p-1)*sg*gp
    rho = F**(-1/p); drho = -(1/p)*F**(-1/p-1)*dF
    return mp.sqrt(rho**2 + drho**2)
def arcs_at(k):
    mp.mp.dps = DPS
    p = 1 + k*H
    a1 = mp.quad(lambda th: integrand(th, p), [0, mp.pi/4])
    a2 = mp.quad(lambda th: integrand(th, p), [mp.pi/4, mp.pi/2])
    a3 = mp.quad(lambda th: integrand(th, p), [mp.pi/2, mp.pi])
    return k, (a1, a2, a3)
if __name__ == '__main__':
    t0 = time.time()
    with Pool(12) as pool:
        vals = dict(pool.map(arcs_at, range(-K, K+1)))
    xs = [k*H for k in range(-K, K+1)]
    M = mp.matrix([[x**j for j in range(2*K+1)] for x in xs])
    res = {}
    for name, sel in [('arc1', lambda v: v[0]), ('total', lambda v: 2*(v[0]+v[1]+v[2]))]:
        ys = mp.matrix([sel(vals[k]) for k in range(-K, K+1)])
        cf = mp.lu_solve(M, ys)
        res[name] = [cf[0], cf[1], 2*cf[2], 6*cf[3]]
    jet = dict(arc1=mp.mpf('5.529732546943809731231005981595168710093'), total=mp.mpf('19.90026636314806698042454547491615647795'))
    # exact arc-1 value from the PSLQ relation: 500 P'''/sqrt5 = 231 z3 + 111 pi^2 log2 + 52 log^3 2 - 80 pi^2 - 1440 log^2 2 + 2400 log2
    z3, l2, pi2 = mp.zeta(3), mp.log(2), mp.pi**2
    ex1 = mp.sqrt(5)/500*(231*z3 + 111*pi2*l2 + 52*l2**3 - 80*pi2 - 1440*l2**2 + 2400*l2)
    print(f"dps={DPS} h={H} K={K} ({time.time()-t0:.0f}s)")
    print(f"  arc1 : FD P'''(1) = {mp.nstr(res['arc1'][3], 45)}")
    print(f"         jet value  = {mp.nstr(jet['arc1'], 45)} ; |FD-jet| = {mp.nstr(abs(res['arc1'][3]-jet['arc1']), 3)}")
    print(f"         PSLQ closed form (231 sqrt5/500 z3 + ...) = {mp.nstr(ex1, 45)} ; |FD-closed| = {mp.nstr(abs(res['arc1'][3]-ex1), 3)}")
    print(f"         same with the EDGE-LAW zeta(3) coefficient (72 sqrt5/500) instead: |FD - that| = {mp.nstr(abs(res['arc1'][3] - (ex1 - mp.sqrt(5)*159/500*z3)), 5)}")
    print(f"         FD P'(1), P''(1) = {mp.nstr(res['arc1'][1], 30)}, {mp.nstr(res['arc1'][2], 30)}")
    print(f"  total: FD P'''(1) = {mp.nstr(res['total'][3], 45)} ; |FD-jet| = {mp.nstr(abs(res['total'][3]-jet['total']), 3)}")
    json.dump({k: [mp.nstr(x, DPS) for x in v] for k, v in res.items()}, open(f'research/10_frontier/E1_edge_law/verify/out/t4_fd_{DPS}_{sys.argv[2]}_{K}.json', 'w'), indent=1)
