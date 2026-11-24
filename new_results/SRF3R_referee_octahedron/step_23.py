import time, json
from mpmath import mp, mpf, quad, sqrt, sin, cos, pi, fabs, log10
t0=time.time()
mp.dps = 55

def S_graph(p):
    """Direct surface quadrature: first-octant graph z=(1-x^p-y^p)^(1/p), raw area element."""
    p = mpf(p)
    def integrand(x, y):
        z = (1 - x**p - y**p)**(1/p)
        zx = -x**(p-1) * z**(1-p)
        zy = -y**(p-1) * z**(1-p)
        return sqrt(1 + zx**2 + zy**2)
    def inner(x):
        ymax = (1 - x**p)**(1/p)
        return ymax * quad(lambda v: integrand(x, ymax*v), [0, 1])
    return 8 * quad(inner, [0, 1])

def S_radial(p):
    """Independent: radial graph over unit sphere, region {x>=y>=0,z>=0} = 1/16 of sphere."""
    p = mpf(p)
    def Nrho(th, ph):
        st, ct = sin(th), cos(th)
        a = (st*cos(ph))**p + (st*sin(ph))**p + ct**p   # = N
        rho = a**(-1/p)
        # d rho: rho = N^{-1/p}
        dNth = p*(st*cos(ph))**(p-1)*ct*cos(ph) + p*(st*sin(ph))**(p-1)*ct*sin(ph) - p*ct**(p-1)*st
        dNph = p*(st*cos(ph))**(p-1)*(-st*sin(ph)) + p*(st*sin(ph))**(p-1)*(st*cos(ph))
        drhoth = -(1/p)*a**(-1/p-1)*dNth
        drhoph = -(1/p)*a**(-1/p-1)*dNph
        st_ = st if st != 0 else mpf(10)**(-mp.dps)  # avoid /0; phi-deriv term ->0 anyway at poles
        grad2 = drhoth**2 + (drhoph/st_)**2
        return rho*sqrt(rho**2 + grad2)*st
    return 16 * quad(lambda ph: quad(lambda th: Nrho(th, ph), [0, pi/2]), [0, pi/4])

def T_simplex(tau):
    tau = mpf(tau)
    f = lambda t1, t2: ((t1*t2*(1-t1-t2))**(-tau)
                        * sqrt(t1**(2*tau) + t2**(2*tau) + (1-t1-t2)**(2*tau)))
    def inner(t1):
        return (1-t1)*quad(lambda v: f(t1, (1-t1)*v), [0,1])
    return quad(inner, [0,1])

def S_via_T(p):
    tau = 1 - 1/mpf(p)
    return 8*(1-tau)**2 * T_simplex(tau)

ref = json.load(open("<project>/research/10_frontier/H1_octahedron_end/direct_S_dps40.json"))
rows=[]
for ds in ["0.01","-0.01","0.03","-0.03"]:
    p = 1+mpf(ds)
    Sg = S_graph(p)
    St = S_via_T(p)
    Sr = S_radial(p)
    h1 = mpf(ref[ds]['N1'])
    d_gt = -log10(fabs(Sg-St)/fabs(St)); d_rt = -log10(fabs(Sr-St)/fabs(St)); d_h1 = -log10(fabs(Sg-h1)/fabs(h1))
    rows.append((ds, Sg, St, Sr, d_gt, d_rt, d_h1))
    print(f"d={ds}: graph-vs-T {float(d_gt):.1f} dig, radial-vs-T {float(d_rt):.1f} dig, graph-vs-H1 {float(d_h1):.1f} dig  [{time.time()-t0:.0f}s]")
for r in rows:
    print(r[0], mp.nstr(r[1], 45))
prog = open("progress.log","a")
for (ds,Sg,St,Sr,dgt,drt,dh1) in rows:
    prog.write(f"[SRF3R] direct quadrature d={ds}: S_graph={mp.nstr(Sg,30)}; graph vs 8(1-tau)^2 T: {float(dgt):.1f} digits; radial-graph vs T: {float(drt):.1f} digits; vs H1: {float(dh1):.1f}\n")
prog.close()
