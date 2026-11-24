# Referee item 5: a_k = [tau^k] S/S(1), k <= K, by an independent exact route.
# f = e^{-tau sum L_i} (3 + X)^{1/2},  X = sum_i (e^{2 tau L_i} - 1),  L_i = log t_i.
# (3+X)^{1/2} = sum_m binom(1/2,m) 3^{1/2-m} X^m  (formal in tau, X = O(tau)),
# X^m = (y1+y2+y3-3)^m with y_i = t_i^{2 tau}, and  int_Delta prod t_i^{a_i} dt1 dt2
#     = prod Gamma(1+a_i) / Gamma(3 + sum a_i)   with a_i = tau(2 j_i - 1).
# So T(tau) = sum_{j} W_j G_j(tau) + O(tau^{K+1}),  G_j = prod Gamma(1+c_i tau)/Gamma(3+c tau),
# c_i = 2 j_i - 1, c = 2J - 3, and log Gamma(1+x) = -gamma x + sum_{n>=2} (-1)^n zeta(n) x^n / n.
import sys, json, itertools
import mpmath as mp

def coeffs(K, dps):
    mp.mp.dps = dps
    lg = [mp.mpf(0), -mp.euler] + [(-1)**n*mp.zeta(n)/n for n in range(2, K+1)]
    fact = [mp.factorial(n) for n in range(K+1)]
    bh = [mp.binomial(mp.mpf(1)/2, m)*mp.power(3, mp.mpf(1)/2-m) for m in range(K+1)]
    T = [mp.mpf(0)]*(K+1)
    for j1 in range(K+1):
        for j2 in range(j1, K+1-j1):
            for j3 in range(j2, K+1-j1-j2):
                J = j1+j2+j3
                mult = len(set(itertools.permutations((j1, j2, j3))))
                W = sum(bh[m]*mp.power(-3, m-J)*fact[m]/(fact[j1]*fact[j2]*fact[j3]*fact[m-J]) for m in range(J, K+1))
                cs = [2*j1-1, 2*j2-1, 2*j3-1]; c = 2*J-3
                A = [mp.mpf(0)] + [lg[n]*(sum(ci**n for ci in cs) - c**n) for n in range(1, K+1)]
                E = [mp.mpf(1)] + [mp.mpf(0)]*K
                for n in range(1, K+1):
                    E[n] = sum(k*A[k]*E[n-k] for k in range(1, n+1))/n
                # divide by (1+c tau)(2+c tau)
                q = [mp.mpf(0)]*(K+1)  # 1/((1+cx)(2+cx)) = sum ((-c)^n (1 - 2^{-n-1})) x^n
                for n in range(K+1):
                    q[n] = mp.power(-c, n)*(1 - mp.power(2, -n-1))
                Wm = W*mult
                for n in range(K+1):
                    T[n] += Wm*sum(E[i]*q[n-i] for i in range(n+1))
    # S/S(1) = (2/sqrt3)(1-tau)^2 T
    a = [2/mp.sqrt(3)*(T[k] - (2*T[k-1] if k >= 1 else 0) + (T[k-2] if k >= 2 else 0)) for k in range(K+1)]
    return T, a

if __name__ == "__main__":
    K = int(sys.argv[1]); dps = int(sys.argv[2])
    T, a = coeffs(K, dps)
    print("T0 - sqrt3/2 =", mp.nstr(T[0]-mp.sqrt(3)/2, 5), "  a0 =", mp.nstr(a[0], 20))
    json.dump({"K": K, "dps": dps, "a": [str(x) for x in a], "T": [str(x) for x in T]},
              open(f"r5_coeffs_K{K}_dps{dps}.json", "w"))
    for k in range(min(K, 8)+1):
        print(k, mp.nstr(a[k], 25))
