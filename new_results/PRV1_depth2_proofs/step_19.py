
import mpmath as mp
mp.mp.dps = 80

def depth2(x, y, p=2):
    rho = max(abs(x), abs(x*y))
    assert 0 < rho < 1
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target*(1-rho))/mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x*y, mp.mpf(0)
    for n in range(1, nmax+1):
        total += h/n**p
        h = x*h + x*xy_power/n
        xy_power *= x*y
    return total

def H_dec(t): return depth2(1-t, -t/(1-t))
def H_dec_closed2(t):
    return (mp.polylog(3, 1-t**2) - mp.polylog(3, 1/(1+t)) + mp.polylog(3, t/(1+t))
            - mp.log(1+t)*mp.polylog(2, -t)/2 - 3*mp.log(1+t)*mp.polylog(2, 1-t)/2
            - mp.log(1+t)*mp.polylog(2, 1-t**2)/4)

print("H_dec closed form #2 (constant C=0):")
for t0 in ['0.1','0.3','0.5','0.61803398874989484820458683436563811772','0.75','0.9']:
    t = mp.mpf(t0)
    print('  t=', t0, ' residual=', mp.nstr(H_dec(t)-H_dec_closed2(t), 8))

phi = (1+mp.sqrt(5))/2; l = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3)
L2 = lambda z: mp.polylog(2, z); L3 = lambda z: mp.polylog(3, z)
print()
print("golden dilogs:")
print("Li2(1/phi^2) - (2z2/5 - l^2) =", mp.nstr(L2(1/phi**2)-(2*z2/5-l**2), 8))
print("Li2(1/phi)   - (3z2/5 - l^2) =", mp.nstr(L2(1/phi)-(3*z2/5-l**2), 8))
print("Li2(-1/phi)  - (-2z2/5 + l^2/2) =", mp.nstr(L2(-1/phi)-(-2*z2/5+l**2/2), 8))
print("golden trilog lattice:")
print("duplication Li3(1/phi^2) - 4Li3(1/phi) - 4Li3(-1/phi) =", mp.nstr(L3(1/phi**2)-4*L3(1/phi)-4*L3(-1/phi), 8))
print("Lewin Li3(1/phi^2) - (4(z3-l z2)/5 + 2l^3/3) =", mp.nstr(L3(1/phi**2)-(4*(z3-l*z2)/5+2*l**3/3), 8))
print("sum rel Li3(1/phi)+Li3(-1/phi) - ((z3-l z2)/5 + l^3/6) =", mp.nstr(L3(1/phi)+L3(-1/phi)-((z3-l*z2)/5+l**3/6), 8))
print()
# specialise closed form 2
t = 1/phi
val = H_dec_closed2(t)
dec1 = 3*z3/5 - 7*l*z2/5 + 2*l**3 + L3(-1/phi) + L3(1/phi)
print("swap_dec closed2 =", mp.nstr(val, 50))
print("DEC1 form        =", mp.nstr(dec1, 50))
print("diff =", mp.nstr(val-dec1, 10))
# my symbolic reduction: Li3(1/phi^2) - 11 l z2/20 + 3 l^3/2
red = L3(1/phi**2) - 11*l*z2/20 + 3*l**3/2
print("closed2 - [Li3(1/phi^2) - 11lz2/20 + 3l^3/2] =", mp.nstr(val - red, 10))
