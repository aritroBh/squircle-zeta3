
import os, json
import mpmath as mp
mp.mp.dps = 60
base = "<project>"
# 2D control: paper law f_n ~ sqrt(2) cos(pi(3-2n)/4) n! (2/pi)^n n^{-3/2} for I_n=[eta^n] I, I=P/(4r(1+eta))
# from H1 n=2 hi coeffs: S_2/S_2(1) in eta=1-1/p-1 = -tau... use hi_coeffs_n2_K80.json
h2 = json.load(open(os.path.join(base,"research/10_frontier/H1_octahedron_end/hi_coeffs_n2_K80.json")))
eta2 = [mp.mpf(x) for x in h2['eta']]  # [eta^k] P/P(1), eta=1-1/p... careful: H1 eta = 1-1/p = tau? For n=2, P/P(1)
# The paper: f_n = I_n with I = P/(4r(1+eta)), eta = 1/p-1 = -tau. So I_n = (-1)^n [tau^n](P/(4sqrt2(1-tau)))
r = mp.sqrt(2)
a = [eta2[k]/(4*r)*(-1)**k for k in range(len(eta2))]  # sign for eta=-tau
# divide by (1+eta)=1/(1-tau)? I = P/(4r(1+eta)) = (P/4r)*(1-eta+...): convolution with (-eta)^k:
I = []
for n in range(len(a)):
    I.append(sum(a[k]*(1)**(n-k) for k in range(n+1)))  # 1/(1+eta)=sum (-eta)^j -> in tau: sum
# law: f_n ~ sqrt2 cos(pi(3-2n)/4) n! (2/pi)^n n^-3/2
print("2D control: ratio I_n / [sqrt2 cos(pi(3-2n)/4) n! (2/pi)^n n^-3/2] should -> 1")
for n in [10,20,40,60,79]:
    law = mp.sqrt(2)*mp.cos(mp.pi*(3-2*n)/4)*mp.factorial(n)*(2/mp.pi)**n*n**(-mp.mpf('1.5'))
    print(f"  n={n}: I_n/law = {mp.nstr(I[n]/law,10)}")
