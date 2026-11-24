# R99: squircle sine near the circle. Own engine: Phi_p(y) = (pi/2) F_p(y)/F_p(1),
# F_p(y) = y 2F1(1/p,1/p;1+1/p;y^p), F_p(1) = pi/(p sin(pi/p)); p-derivatives by mpmath.taylor.
import sys, time
from mpmath import mp, mpf, hyp2f1, sin, cos, pi, log, zeta, gamma, quad, taylor, nstr, diff, floor, fsum
mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
def Phi(p, y):
    return (pi/2)*y*hyp2f1(1/p, 1/p, 1+1/p, y**p)/(pi/(p*sin(pi/p)))
cache = {}
def Psi(phi):
    if phi not in cache:
        y = sin(phi)
        cache[phi] = taylor(lambda p: Phi(p, y), mpf(2), 3)  # [Phi, Phi', Phi''/2, Phi'''/6]
    return cache[phi]
def bder(n):
    # b_n(2+d) = (4/(n pi)) int_0^{pi/2} cos(n Phi_{2+d}(sin f)) cos f df ; Taylor in d up to d^3
    def comp(f, k):
        t = Psi(f); A = n*t[0]; x1, x2, x3 = n*t[1], n*t[2], n*t[3]
        cA, sA = cos(A), sin(A)
        # cos(A+X), X = x1 d + x2 d^2 + x3 d^3
        c = [cA, -sA*x1, -sA*x2 - cA*x1**2/2, -sA*x3 - cA*x1*x2 + sA*x1**3/6]
        return c[k]*cos(f)
    out = []
    for k in range(4):
        v = 4/(n*pi)*quad(lambda f: comp(f, k), [0, pi/8, pi/4, 3*pi/8, pi/2])
        out.append(v*[1, 1, 2, 6][k])
    return out  # b_n, b_n', b_n'', b_n'''
t0 = time.time()
l2 = log(2)
for n in [1, 3, 5, 7, 9, 11]:
    b = bder(n)
    if n == 1:
        pr1 = (l2-1)/2
        pr2 = mpf(3)/8 - l2 + l2**2/4 + 5*pi**2/96 - pi**4/768
        print("n=1: b=%s  b'=%s  [printed (log2-1)/2: diff %s]" % (nstr(b[0], 20), nstr(b[1], 25), nstr(b[1]-pr1, 3)))
        print("     b''=%s  [notes' closed form: diff %s]" % (nstr(b[2], 25), nstr(b[2]-pr2, 3)))
        print("     b'''=%s  [notes: -0.30262285104154317677477849746]" % nstr(b[3], 25))
    else:
        fl = 2*floor(mpf(n-1)/4)+1
        pr1 = mpf((-1)**((n-1)//2))/((n*n-1)*fl)
        print("n=%d: b'=%s  printed law %s  diff %s ; b''=%s ; b'''=%s" % (n, nstr(b[1], 22), nstr(pr1, 10), nstr(b[1]-pr1, 3), nstr(b[2], 22), nstr(b[3], 22)))
    sys.stdout.flush()
# b_3'' printed in notes: (12-16 ln2 + 3 pi^2)/128
print("time", round(time.time()-t0, 1), "s")
