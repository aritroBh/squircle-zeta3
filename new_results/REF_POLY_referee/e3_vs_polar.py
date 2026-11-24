"""Check the referee's exact edge expansion E_k(t) (edge_sym.py) against the polar jet, arc by arc."""
import sympy as sp, pickle, mpmath as mp
from edge_sym import t, L, al, be, B1, C1, x0, y0, x1, y1, X0, X1, Xf1, Xf2
from jet import integrand, kinks
from polys import polys
from geom import edges, edge_data, kink_of
mp.mp.dps = 45
E = pickle.load(open('edge_E_general.pkl', 'rb'))
syms = (t, L, al, be, B1, C1, x0, y0, x1, y1, X0, X1, Xf1, Xf2)
fs = [sp.lambdify(syms, E[k], 'mpmath') for k in range(4)]
def sector_jet(U, i, N=4, h=mp.mpf(1)/128):
    ks = kinks(U); m = len(ks)
    tha, ja = ks[i % m]; thb, jb = ks[(i+1) % m]
    if thb <= tha: thb += mp.pi
    Wd = thb-tha; acc = [mp.mpf(0)]*N; k = 0; cut = mp.mpf(10)**(-(mp.mp.dps+5)); W = [1]*len(U)
    while True:
        done = True
        for tt in ([k*h] if k == 0 else [k*h, -k*h]):
            q = mp.exp(-mp.pi*mp.sinh(tt)); s = 1/(1+q); sc = q/(1+q)
            if min(s, sc) < cut: continue
            done = False; wgt = mp.pi*mp.cosh(tt)*s*sc*Wd*h
            da = Wd*s; db = Wd*sc; th = tha+da if s < sc else thb-db
            f = integrand(th, da, db, ja, jb, U, W, N)
            for n in range(N): acc[n] += wgt*f[n]
        if done and k > 0: break
        k += 1
    return acc, tha, thb
for name in ['gen_hexagon', 'gen_octagon', 'reg_octagon']:
    U = polys()[name]
    ks = kinks(U)
    for i in range(2):
        acc, tha, thb = sector_jet(U, i)
        # the edge whose corners lie on the rays tha, thb
        r = lambda th, j: (lambda d_: ((-U[j][1]*mp.cos(0)) , 0))
        V = []
        for th, j in ((tha, ks[i][1]), (thb, ks[(i+1) % len(ks)][1])):
            dvec = (mp.cos(th), mp.sin(th))
            nrm = mp.fsum(abs(U[q][0]*dvec[0]+U[q][1]*dvec[1]) for q in range(len(U)) if q != j)
            V.append((dvec[0]/nrm, dvec[1]/nrm))
        A, B = V
        Le, h, L0, L1, xi0, xi1, j0, j1 = edge_data(U, A, B)
        far = [q for q in range(len(U)) if q not in (j0, j1)]
        aff = lambda q: (lambda tt: abs(U[q][0]*(A[0]+tt*(B[0]-A[0]))+U[q][1]*(A[1]+tt*(B[1]-A[1]))))
        a1 = aff(far[0]); b1 = a1(0); c1 = a1(1)-a1(0)
        a2 = aff(far[1]) if len(far) > 1 else None
        def Ek(k, tt):
            Xf2v = mp.log(a2(tt)) if a2 else 0
            return fs[k](tt, Le, xi0, xi1, b1, c1, A[0], A[1], B[0], B[1],
                         mp.log(xi0*tt), mp.log(xi1*(1-tt)), mp.log(a1(tt)), Xf2v)
        res = [mp.quad(lambda tt: Ek(k, tt), [0, mp.mpf(1)/2, 1]) for k in range(4)]
        print(name, 'arc', i, ' |int E_k dt - polar sector coefficient|, k=0..3:',
              ['%.1e' % float(abs(res[k]-acc[k])) for k in range(4)], flush=True)
