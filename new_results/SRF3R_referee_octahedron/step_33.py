from mpmath import mp, mpf, quad, sqrt, exp, pi, log
mp.dps = 40
z0 = (-log(2)+pi*1j)/2
R = abs(z0)

def Phi(tau, V):
    tau=mpf(tau); V=mpf(V)
    em = exp(-V)
    return (exp(tau*V)*(1-em)**(1-2*tau)
            * quad(lambda s: ((s*(1-s))**(-tau)
                   * sqrt((1-em)**(2*tau)*(s**(2*tau)+(1-s)**(2*tau)) + exp(-2*tau*V))), [0,1]))

print("w fixed, V large: compare Phi(w/V, V) with candidates")
for w, V in [(0.3, 30), (0.3, 100), (0.7, 100)]:
    tau = mpf(w)/V
    ph = Phi(tau, V)
    H_srf3 = exp(w)*sqrt(1+2*exp(2*w))      # proof.md's H0 = e^w sqrt(1+2e^{2w})
    H_true = sqrt(1+2*exp(2*w))            # my candidate = e^w sqrt(2+e^{-2w})
    H_alt  = exp(w)*sqrt(2+exp(-2*w))
    print(f"w={w}, V={V}: Phi={mp.nstr(ph,12)}")
    print(f"   e^w sqrt(1+2e^(2w)) = {mp.nstr(H_srf3,12)}  ratio {mp.nstr(ph/H_srf3,8)}")
    print(f"   sqrt(1+2e^(2w))     = {mp.nstr(H_true,12)}  ratio {mp.nstr(ph/H_true,8)}")
    print(f"   e^w sqrt(2+e^(-2w)) = {mp.nstr(H_alt,12)}  ratio {mp.nstr(ph/H_alt,8)}")

# strength mismatch check: per-edge strength from S4.1 vs strength of sqrt(1+2e^{2z}) at z0
# H(z) = sqrt(1+2e^{2z}): H^2 ~ 4 e^{2 z0} (z-z0) = -2 (z-z0); so H ~ sqrt(2) (z0-z)^{1/2} up to branch
# i.e. H ~ i sqrt(2) z0^{1/2} (1-z/z0)^{1/2} (one determination)
import mpmath as m
str_H = 1j*sqrt(2)*z0**mpf('0.5')
str_edge = (pi*1j/4)*exp(z0)*z0**mpf('0.5')
print("\nstrength of sqrt(1+2e^{2z}) at z0 (one branch):", mp.nstr(str_H,10))
print("per-edge strength from S4.1 (pi i/4) e^{z0} z0^{1/2}:", mp.nstr(str_edge,10))
print("ratio:", mp.nstr(str_edge/str_H,8))
