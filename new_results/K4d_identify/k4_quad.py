
"""tanh-sinh quadrature for 1D and nested 2D integrals."""
import mpmath as mp
from mpmath import mpf

def tanh_sinh_nodes(a, b, h, T, wmin=None):
    """nodes (x,w) for int_a^b f dx; t = k h, |t|<=T.
    K4b fix: nodes with |w| <= wmin are skipped (their total mass is
    negligible and their x-coordinates can carry no valid digits)."""
    nodes = []
    k = 0
    while True:
        t = k*h
        if t > T: break
        sh = mp.sinh(t)
        u = mp.tanh(mp.pi/2*sh)
        x = (a+b)/2 + (b-a)/2*u
        # dx/dt
        ch = mp.cosh(t)
        cp = mp.cosh(mp.pi/2*sh)
        w = (b-a)/2 * (mp.pi/2*ch)/(cp*cp) * h
        nodes.append((x, w, k))
        k += 1
    # symmetrize: k and -k give symmetric nodes
    out = []
    for (x,w,k) in nodes:
        if wmin is not None and abs(w) <= wmin:
            continue
        if k == 0:
            out.append((x,w))
        else:
            if x > a and x < b:
                out.append((x,w))
            x2 = a+b-x
            if x2 > a and x2 < b:
                out.append((x2, w))
    return out

def ts_quad(f, a, b, h, T):
    s = mpf(0)
    for (x,w) in tanh_sinh_nodes(a,b,h,T):
        fx = f(x)
        if fx != 0:
            s += w*fx
    return s
