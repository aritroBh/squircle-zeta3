import mpmath as mp, pickle
from fractions import Fraction as Fr
mp.mp.dps = 100
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator

def Li12(x, y, N=2000):
    # standard Li_{1,2}(x,y) = sum_{0<m<n} x^m y^n/(m n^2)
    x = M(x); y = M(y); s = mp.mpf(0)
    yn = y
    for n in range(2, N):
        yn *= y
        t = mp.mpf(0); xm = x
        for m in range(1, n):
            t += xm/m; xm *= x
        s += t*yn/n**2
        if abs(yn) < mp.mpf(10)**(-mp.mp.dps-10): break
    return s

Lam = {k: mp.mpf(v) for k, v in pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"].items()}
tests = {"Lam1": (Fr(1,3), Fr(-1,3)), "Lam2": (Fr(3,1), Fr(-1,9)),
         "Lam3": (Fr(2,5), Fr(1,4)),  "Lam4": (Fr(5,2), Fr(1,10))}
for name,(x,y) in tests.items():
    v = Li12(x,y)
    print(name, "Li_{1,2}(%s,%s) ="%(x,y), mp.nstr(v,25))
    print("   stored:", mp.nstr(Lam[name],25), " diff:", mp.nstr(abs(v-Lam[name]),3))

# stuffle numeric check: Li2(x)Li1(y) = Li21A(x,y) + Li12(y,x) + Li3(xy)
x, y = 0.5, 0.5
def Li21A(x, y, N=800):
    x = mp.mpf(x); y = mp.mpf(y); s = mp.mpf(0); yn = y
    for n in range(2, N):
        yn *= y
        t = mp.mpf(0); xm = x
        for m in range(1, n):
            t += xm/m**2; xm *= x
        s += t*yn/n
    return s
lhs = mp.polylog(2,x)*(-mp.log(1-y))
rhs = Li21A(x,y) + Li12(y,x) + mp.polylog(3, x*y)
print("stuffle: lhs-rhs =", mp.nstr(abs(lhs-rhs),3))
