# Planner check of HEX5: regular hexagon K_p = {sum_{i=1..3} |<n_i,x>|^p <= 1}, n_i at angles i*pi/3; P^(5)(1) by polar FD + Richardson,
# and HEX5's closed form 128 P5/sqrt3 = ... evaluated directly.
import mpmath as mp
mp.mp.dps = 70
N = [(mp.cos(i*mp.pi/3), mp.sin(i*mp.pi/3)) for i in (1,2,3)]
def P(p):
    def g(th):
        s, c = mp.sin(th), mp.cos(th); F = 0; Fp = 0
        for (a,b) in N:
            v = a*c + b*s; dv = -a*s + b*c
            F += abs(v)**p; Fp += p*abs(v)**(p-1)*mp.sign(v)*dv
        rho = F**(-1/p); drho = -(1/p)*F**(-1/p-1)*Fp
        return mp.sqrt(rho**2+drho**2)
    kinks = [mp.pi/2 + i*mp.pi/3 for i in range(-1,3)]   # directions where some <n_i,e_theta> = 0, within [0, pi]
    br = sorted(set([mp.mpf(0)] + [k % mp.pi for k in kinks] + [mp.pi]))
    return 2*sum(mp.quad(g,[br[i],br[i+1]]) for i in range(len(br)-1))
print("P(1) =", mp.nstr(P(mp.mpf(1)),20), " 2sqrt3 =", mp.nstr(2*mp.sqrt(3),20))
def d5(h):
    f = {j: P(1+j*h) for j in (-3,-2,-1,1,2,3)}
    return (-f[-3]+4*f[-2]-5*f[-1]+5*f[1]-4*f[2]+f[3])/(2*h**5)
D = [d5(h) for h in (mp.mpf('4e-4'), mp.mpf('2e-4'), mp.mpf('1e-4'))]
R1 = [(4*D[i+1]-D[i])/3 for i in range(2)]; R2 = (16*R1[1]-R1[0])/15
l = mp.log(2); pi = mp.pi; z3 = mp.zeta(3); z5 = mp.zeta(5); z2 = pi**2/6
closed = mp.sqrt(3)/128*(30720*l-4920*pi**2-61440*l**2+6240*pi**2*l+30720*l**3+448*pi**4+5760*z3*l-1920*pi**2*l**2-5120*l**4
         +11340*z5+5580*z2*z3-123*pi**4*l-1440*z3*l**2+160*pi**2*l**3+256*l**5)
print("FD Richardson P5 =", mp.nstr(R2,18)); print("HEX5 closed form =", mp.nstr(closed,18)); print("HEX5 value 656.41536196762475447664760905831 ; diff FD-closed =", mp.nstr(R2-closed,3))
