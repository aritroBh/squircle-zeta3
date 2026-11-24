# Planner check of HEXC2: (1) count, by own exact enumeration, the Mobius maps sending the five points {0, 1/2, 1, -1, inf}
# injectively into the sixth roots of unity; (2) the task map's images.
import itertools, sympy as sp
z = sp.exp(sp.I * sp.pi / 3)
roots = [sp.nsimplify(sp.expand_complex(z**j)) for j in range(6)]
pts = [0, sp.Rational(1, 2), 1, -1, sp.oo]
def mobius_from(p, q):  # map sending 0, 1, oo to p, q, r is determined by 3 points; enumerate images of (0, 1, oo), then check 1/2 and -1
    pass
count = 0; maps = []
for a0, a1, ainf in itertools.permutations(range(6), 3):
    p, q, r = roots[a0], roots[a1], roots[ainf]
    # f(z) = (r*(z)*(q - p) ... ) : the Mobius map with f(0)=p, f(1)=q, f(oo)=r is f(z) = (r*z*(q-p) + p*(r-q)) / (z*(q-p) + (r-q))
    f = lambda x: sp.simplify((r * x * (q - p) + p * (r - q)) / (x * (q - p) + (r - q)))
    ih, im1 = f(sp.Rational(1, 2)), f(-1)
    def idx(v):
        for j, w in enumerate(roots):
            if sp.simplify(sp.expand_complex(v - w)) == 0: return j
        return None
    jh, jm = idx(ih), idx(im1)
    if jh is not None and jm is not None and len({a0, a1, ainf, jh, jm}) == 5:
        count += 1; maps.append((a0, jh, a1, jm, ainf))
print("Mobius maps sending {0, 1/2, 1, -1, inf} injectively into mu_6:", count)
print("images (0, 1/2, 1, -1, inf) as powers of zeta:", maps)
f = lambda x: (z**2 * x + 1) / (z**4 * x + 1)
print("task map f(z) = (z^2 x + 1)/(z^4 x + 1): f(0), f(1/2), f(1), f(-1) =", [sp.nsimplify(sp.simplify(sp.expand_complex(f(x) / z**k))) for x, k in [(0, 0), (sp.Rational(1, 2), 1), (1, 2), (-1, 5)]], "(each should be 1); f(oo) = z^2/z^4 = z^-2 = z^4")
