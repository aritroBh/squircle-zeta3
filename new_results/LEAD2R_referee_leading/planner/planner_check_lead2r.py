# Planner check: the leading-coefficient law c(w) = eps_{w-2} A_{w-1} [y^w]E(y) / (2 (w-2)!), E(y) = Gamma(1+y/2)Gamma(1/2-y/2)/sqrt(pi),
# against independently known values: c(3) = -(z3/8 + pi^2 l/24 + l^3/12) (weight-3 law), c(4) = D2/2, c(5) = E3/6, c(6) from the kernels.
import mpmath as mp
mp.mp.dps = 50
E = lambda y: mp.gamma(1 + y/2)*mp.gamma(mp.mpf(1)/2 - y/2)/mp.sqrt(mp.pi)
co = mp.taylor(E, 0, 12)
def zigzag(n):  # Euler zigzag numbers A_n via the Seidel-Entringer triangle
    row = [1]
    for i in range(1, n + 1):
        new = [0]
        for x in reversed(row): new.append(new[-1] + x)
        row = new
    return row[-1]
eps = lambda m: [1, -1, -1, 1][m % 4]
c = lambda w: eps(w - 2)*zigzag(w - 1)*co[w]/(2*mp.factorial(w - 2))
l = mp.log(2); z3 = mp.zeta(3); pi = mp.pi
print("zigzag A_0..A_7:", [zigzag(n) for n in range(8)])
print("c(3) - [-(z3/8 + pi^2 l/24 + l^3/12)] =", mp.nstr(c(3) + (z3/8 + pi**2*l/24 + l**3/12), 5))
D2 = -l**4/24 - pi**2*l**2/24 - pi**4/160 - z3*l/4
print("c(4) - D2/2 =", mp.nstr(c(4) - D2/2, 5))
print("c(5) - 0.409481194212564356915219215067 =", mp.nstr(c(5) - mp.mpf('0.409481194212564356915219215067'), 5))
print("c(6) - 0.335687558224020357007672279844 =", mp.nstr(c(6) - mp.mpf('0.335687558224020357007672279844'), 5))
print("c(2) =", mp.nstr(c(2), 20))
