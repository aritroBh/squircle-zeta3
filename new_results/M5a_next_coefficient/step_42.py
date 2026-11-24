import mpmath as mp
mp.mp.dps = 80
r2 = mp.sqrt(2); b = r2 - 1; s2i = 1/r2
def li_mn(m, n, x, y):
    # series sum_{p>q>=1} x^p y^q/(p^m q^n); needs |x*y|<=1, |x|<1-ish; truncate adaptively
    tot = mp.mpf(0)
    p = 1
    while True:
        row = mp.fsum(x**p * y**q / (p**m * q**n) for q in range(1, p))
        tot += row
        if abs(row) < mp.mpf('1e-85') and p > 50:
            break
        p += 1
        if p > 3000: break
    return tot
A33_b1  = li_mn(3, 3, b, 1)
A33_s1  = li_mn(3, 3, s2i, 1)
A33_b2  = li_mn(3, 3, b, 2)
A42_bs  = li_mn(4, 2, b, s2i)
print("Li_{3,3}(b,1)     =", mp.nstr(A33_b1, 60))
print("Li_{3,3}(1/r2,1)  =", mp.nstr(A33_s1, 60))
print("Li_{3,3}(b,2)     =", mp.nstr(A33_b2, 60))
print("Li_{4,2}(b,1/r2)  =", mp.nstr(A42_bs, 60))
N6 = (-12 + 209*r2)*A33_b1 + (-160 + 80*r2)*A33_s1 + (mp.mpf(329)/4*r2)*A33_b2 + 120*r2*A42_bs
print("N6 (numeric)      =", mp.nstr(N6, 60))
W6 = mp.mpf(open('<project>/research/08_theory/M5a_next_coefficient/W6_400.txt').read().strip())
res = W6 - N6
print("W6 - N6           =", mp.nstr(res, 60))
import json
json.dump({"Li33_b1": mp.nstr(A33_b1, 75), "Li33_s1": mp.nstr(A33_s1, 75),
           "Li33_b2": mp.nstr(A33_b2, 75), "Li42_bs": mp.nstr(A42_bs, 75),
           "N6": mp.nstr(N6, 75), "W6_minus_N6": mp.nstr(res, 75)},
          open('<project>/research/08_theory/M5a_next_coefficient/n6_values.json', 'w'), indent=1)
print("saved n6_values.json")