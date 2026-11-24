
import mpmath as mp
mp.mp.dps = 50

def phi(s, th):
    c = mp.cos(th)
    return 1 - (s + c)/mp.sqrt(1 + 2*c*s + s**2)

def S_of_U(u, th):
    C = mp.cos(th/2)**2; Z = mp.sin(th/2)**2
    return C/u - mp.cos(th) - Z*u

def Lambda(nu, th):
    z = mp.sin(th/2)**2
    return 2*z*mp.hyper([-nu, 1+nu], [2], z)

print("Corrected Mellin: E[S^w] = (pi w/sin pi w) Lambda_w/(1-c)")
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3, mp.pi/6]:
    c = mp.cos(th)
    for w in [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.77')]:
        lhs = mp.quad(lambda s: s**w*phi(s,th), [0, mp.inf])/(1-c)
        lhs2 = mp.quad(lambda u: S_of_U(u, th)**w, [0,1])
        rhs = (mp.pi*w/mp.sin(mp.pi*w))*Lambda(w, th)/(1-c)
        d = max(abs(lhs-rhs), abs(lhs2-rhs))
        print(f" th={mp.nstr(th,6)} w={w}: quad={mp.nstr(lhs,13)} explaw={mp.nstr(lhs2,13)} leg={mp.nstr(rhs,13)} maxdiff={mp.nstr(d,2)}")

# Cumulants of log(1/S_theta) from the law, vs (-1)^j j! a_j from the paper's expansion
# log R_theta(eps) = log(2z) + sum a_j eps^j ; R = (pi eps/sin pi eps) Lambda_eps
def a_coeffs(th, J=4):
    z = mp.sin(th/2)**2
    f = lambda e: mp.log((mp.pi*e/mp.sin(mp.pi*e))*Lambda(e, th)/(2*z))
    return [mp.diff(f, mp.mpf(0), j)/mp.factorial(j) for j in range(1, J+1)]

def cumulants_numeric(th, J=4):
    # K(t) = log E[exp(-t log S)] = log E[S^{-t}]; cumulants of L=log(1/S): k_j = K^{(j)}(0)
    f = lambda t: mp.log(mp.quad(lambda u: S_of_U(u, th)**(-t), [0,1]))
    return [mp.diff(f, mp.mpf(0), j) for j in range(1, J+1)]

print("\nCumulants k_j vs (-1)^j j! a_j:")
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    aa = a_coeffs(th)
    kk = cumulants_numeric(th)
    for j in range(4):
        target = (-1)**(j+1)*mp.factorial(j+1)*aa[j]
        print(f" th={mp.nstr(th,6)} k{j+1}: numeric={mp.nstr(kk[j],14)}  closed={mp.nstr(target,14)}  diff={mp.nstr(abs(kk[j]-target),2)}")
