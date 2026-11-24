
import mpmath as mp
mp.mp.dps = 50

# ---------- Task 1: the theta-dependent slope law ----------
# phi_theta(s) = d/ds h_c(s), h_c(s) = 1 + s - sqrt(1 + 2 c s + s^2), c = cos(theta)
def phi(s, th):
    c = mp.cos(th)
    return 1 - (s + c)/mp.sqrt(1 + 2*c*s + s**2)

def hc(s, th):
    c = mp.cos(th)
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)

for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    c = mp.cos(th)
    mass = mp.quad(lambda s: phi(s, th), [0, mp.inf])
    print(f"theta={mp.nstr(th,8)}  int phi = {mp.nstr(mass,20)}   1-cos = {mp.nstr(1-c,20)}  diff={mp.nstr(abs(mass-(1-c)),3)}")

# Geometric identity: phi_theta(s) = 1 - cos(beta(s)), beta = angle between chord u+sv and v
# and the change of variables s = sin(theta-beta)/sin(beta), tan(beta/2) = tan(theta/2)*U
# => S_theta = cos^2(theta/2)/U - cos(theta) - sin^2(theta/2)*U, U uniform(0,1)
def S_of_U(u, th):
    C = mp.cos(th/2)**2; Z = mp.sin(th/2)**2
    return C/u - mp.cos(th) - Z*u

print()
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    # verify density of S_of_U equals phi_theta: P(S>s) = measure of {u: S(u)>s}; compare with tail
    for s0 in [mp.mpf('0.3'), 1, 3]:
        # tail from density
        tail = mp.quad(lambda s: phi(s,th), [s0, mp.inf])
        # U(s): solve S(u)=s -> C/u - Z u - c = s -> Z u^2 + (s+c) u - C = 0
        Z = mp.sin(th/2)**2; C = mp.cos(th/2)**2; c = mp.cos(th)
        u_star = (-(s0+c) + mp.sqrt((s0+c)**2 + 4*Z*C))/(2*Z)
        print(f"th={mp.nstr(th,6)} s={s0}: tail={mp.nstr(tail,15)}  u*={mp.nstr(u_star,15)}  diff={mp.nstr(abs(tail-u_star),3)}")
