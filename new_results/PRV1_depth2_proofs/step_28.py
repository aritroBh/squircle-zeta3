
import mpmath as mp
mp.mp.dps = 60
phi = (1+mp.sqrt(5))/2
def depth2(x, y, p=2, q=1):
    # Li_{p,q}(x,y) = sum x^n y^m / (n^p m^q)
    rho = max(abs(x), abs(x*y)); assert 0 < rho < 1
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target*(1-rho))/mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x*y, mp.mpf(0)
    for n in range(1, nmax+1):
        total += h/n**p
        h = x*h + x*xy_power/n**(2-q+1) if False else x*h + x*xy_power/n
        xy_power *= x*y
    return total
# simpler: direct generic
def Li_pq(x, y, p, q):
    rho = max(abs(x), abs(x*y)); assert 0 < rho < 1
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target*(1-rho))/mp.log(rho)))
    total = mp.mpf(0)
    xn = x
    for n in range(1, nmax+2):
        # inner sum over m=1..n-1 of y^m/m^q
        s = mp.mpf(0); ym = y
        for m in range(1, n):
            s += ym/m**q; ym *= y
        total += xn/n**p * s
        xn *= x
    return total

v = Li_pq(-1/phi, -1/phi, 1, 2)
print("Li_{1,2}(-1/phi,-1/phi) =", mp.nstr(v, 50))
print("stored Li12_11          = -0.09051350756610804496858303299057246654601003439637")
# stuffle check: Lam + Li12 = Li2(x)Li1(y) - Li3(xy)
x, y = -1/phi, -1/phi
lam = Li_pq(x, y, 2, 1)
stuffle = mp.polylog(2,x)*mp.polylog(1,y) - mp.polylog(3, x*y)
print("Lam + Li12 - (Li2Li1 - Li3) =", mp.nstr(lam + v - stuffle, 10))
# my reduced Lambda vs series
l = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3)
Lam_red = z3 + l*z2/5 - l**3/2 - 2*mp.polylog(3, 1/phi)
print("Lam_red - series =", mp.nstr(Lam_red - lam, 10))
swap_red = 4*z3/5 - 8*l*z2/5 + 13*l**3/6
sw = Li_pq(1/phi**2, -phi, 2, 1)
print("swap_red - series =", mp.nstr(swap_red - sw, 10))
# DEC1's forms vs series
Lam_DEC1 = 2*z3/5 + 4*l*z2/5 - l**3 + 3*mp.polylog(3,-1/phi) + mp.polylog(3,1/phi)
swap_DEC1 = 3*z3/5 - 7*l*z2/5 + 2*l**3 + mp.polylog(3,-1/phi) + mp.polylog(3,1/phi)
print("Lam_DEC1 - series =", mp.nstr(Lam_DEC1 - lam, 10))
print("swap_DEC1 - series =", mp.nstr(swap_DEC1 - sw, 10))
