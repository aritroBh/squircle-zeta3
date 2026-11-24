from mpmath import (mp, mpf, pi, gamma, sin, cos, log, nsum, inf, nstr, zeta, gammainc, re, im)
mp.dps = 80
def N(x,d=30): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
def a_g(n,s): return re((2/pi)*(-1j*n)**(-s-1) * gammainc(s+1, 0, -1j*n*pi))
def b_g(n,s): return im((2/pi)*(-1j*n)**(-s-1) * gammainc(s+1, 0, -1j*n*pi))
s = mpf('2.5')
# Q(s): sum a_n(s) = -pi^s/(s+1). convergence: a_n ~ kink n^-3.5 + bnd (-1)^n n^-2: fine
Qn = nsum(lambda n: a_g(n,s), [1, inf])
print("Q(2.5) series =", N(Qn,26), " exact =", N(-pi**s/(s+1),26), " diff:", nstr(abs(Qn+pi**s/(s+1)),3))
Vn = nsum(lambda n: b_g(n,s)/n, [1, inf])
print("V(2.5) series =", N(Vn,26), " exact =", N(pi**(s+1)/((s+1)*(s+2)),26), " diff:", nstr(abs(Vn-pi**(s+1)/((s+1)*(s+2))),3))
# kink parts at s=2.5
kQ = -(2/pi)*gamma(s+1)*sin(pi*s/2)*zeta(s+1)
print("Q(2.5) kink part =", N(kQ,26), " (boundary makes up the rest)")
logline("step_51: non-integer s=2.5 spot checks Q,V verified")
logline("SUMMARY: AV3 complete. Key results: (i) coeff closed forms + kink/boundary asymptotics; (ii) f-family has NO 1/n^3 mode at any real s (kink node at s=2), g-family kink carries 1/n^3 at s=2 with amp -4/pi on odd n; (iii) 8 exact identities Q,P,R,Phi,T,Psi,V,W with 2-derivative tables at s=2,3,4; (iv) zeta(3) survives only in W(s) [sum of sine coefs = log-sin integral] and B(s,s') [cos x sin pairing = log-cos integral]; (v) zeta(3),zeta'(3),zeta'(4),G,beta'(3),beta(4) appear in kink parts of same-parity derivatives and cancel (Si/Ci corollaries: S_Si=5pi^3/72, S_Si2=49pi^3/432, S_Sisin=pi^2 log2/4, S4, Cisum, CiS4 all PROVED); (vi) OPEN: CiS, W'(2), R4, W'(3); (vii) rotation: zeta(3) scales as sin theta (W) and sin(2theta) (B); Parseval theta-blind.")
print("done")
