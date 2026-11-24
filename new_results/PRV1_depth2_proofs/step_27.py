
import pickle, mpmath as mp
d = pickle.load(open('research/08_theory/DEC1_decagon/li12_830.pkl','rb'))
mp.mp.dps = 60
phi = (1+mp.sqrt(5))/2
for k, v in d.items():
    print(k, '=', mp.nstr(v, 50))
# Li12_11 should be Lambda_dec = Li_21(-1/phi,-1/phi); Li12_mphi_x4 = swap = Li_21(1/phi^2, -phi)
def depth2(x, y, p=2):
    rho = max(abs(x), abs(x*y)); assert 0 < rho < 1
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target*(1-rho))/mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x*y, mp.mpf(0)
    for n in range(1, nmax+1):
        total += h/n**p
        h = x*h + x*xy_power/n
        xy_power *= x*y
    return total
lam = depth2(-1/phi, -1/phi); sw = depth2(1/phi**2, -phi)
print("series Lam =", mp.nstr(lam, 50))
print("series swap=", mp.nstr(sw, 50))
print("diff Lam:", mp.nstr(d['Li12_11']-lam, 10))
print("diff swap:", mp.nstr(d['Li12_mphi_x4']-sw, 10))
# now the fully reduced forms
l = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3)
Lam_red = z3 + l*z2/5 - l**3/2 - 2*mp.polylog(3, 1/phi)
swap_red = 4*z3/5 - 8*l*z2/5 + 13*l**3/6
print("Lam reduced - 830 value =", mp.nstr(Lam_red - d['Li12_11'], 10))
print("swap reduced - 830 value =", mp.nstr(swap_red - d['Li12_mphi_x4'], 10))
print("swap reduced =", mp.nstr(swap_red, 50))
