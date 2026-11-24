
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 150
import rul1, rul1num

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

l2 = mp.log(2)
z = {n: mp.zeta(n) for n in range(2, 10)}
Li = {n: mp.polylog(n, mp.mpf(1)/2) for n in range(4, 10)}
ex = pickle.load(open('extras_val.pkl','rb'))
X6, X7a, X7b = ex['zbar51'], ex['zbar511'], ex['zbar421']

def I_word(m, j, dps=150):
    w = (0,)*(m+1) + (-1,)*j
    return mp.factorial(j)*mp.factorial(m) * rul1num.g_word_over_0m1_num(w, dps=dps)

basis_vals = {
 2: [('z2', z[2]), ('log2^2', l2**2)],
 3: [('z3', z[3]), ('z2*log2', z[2]*l2), ('log2^3', l2**3)],
 4: [('z2^2', z[2]**2), ('z3*log2', z[3]*l2), ('z2*log2^2', z[2]*l2**2), ('log2^4', l2**4), ('Li4', Li[4])],
 5: [('z5', z[5]), ('z2*z3', z[2]*z[3]), ('z3*log2^2', z[3]*l2**2), ('z2^2*log2', z[2]**2*l2), ('z2*log2^3', z[2]*l2**3),
     ('log2^5', l2**5), ('Li5', Li[5]), ('Li4*log2', Li[4]*l2)],
 6: [('log2^6', l2**6), ('z2*log2^4', z[2]*l2**4), ('z2^2*log2^2', z[2]**2*l2**2), ('z2^3', z[2]**3),
     ('z3*log2^3', z[3]*l2**3), ('z3*z2*log2', z[3]*z[2]*l2), ('z3^2', z[3]**2), ('z5*log2', z[5]*l2),
     ('Li4*log2^2', Li[4]*l2**2), ('Li4*z2', Li[4]*z[2]), ('Li5*log2', Li[5]*l2), ('Li6', Li[6]), ('X6', X6)],
 7: [('log2^7', l2**7), ('z2*log2^5', z[2]*l2**5), ('z2^2*log2^3', z[2]**2*l2**3), ('z2^3*log2', z[2]**3*l2),
     ('z3*log2^4', z[3]*l2**4), ('z3*z2*log2^2', z[3]*z[2]*l2**2), ('z3*z2^2', z[3]*z[2]**2), ('z3^2*log2', z[3]**2*l2),
     ('z5*log2^2', z[5]*l2**2), ('z5*z2', z[5]*z[2]), ('Li4*log2^3', Li[4]*l2**3), ('Li4*z2*log2', Li[4]*z[2]*l2),
     ('Li4*z3', Li[4]*z[3]), ('Li5*log2^2', Li[5]*l2**2), ('Li5*z2', Li[5]*z[2]), ('Li6*log2', Li[6]*l2),
     ('X6*log2', X6*l2), ('z7', z[7]), ('Li7', Li[7]), ('X7a', X7a), ('X7b', X7b)],
}

reds = {}
for n in range(2, 8):
    bvals = basis_vals[n]
    for m in range(n-1):
        j = n-1-m
        Iv = I_word(m, j)
        vec = [Iv] + [v for _, v in bvals]
        rel = mp.pslq(vec, tol=mp.mpf('1e-130'), maxcoeff=10**14, maxsteps=10000)
        if rel is None:
            print(f"I_{{{m},{j}}} (w{n}): PSLQ FAIL")
            continue
        a0 = rel[0]
        coeffs = [Fr(-x, a0) for x in rel[1:]]
        # verify exactly the numeric relation
        res = Iv + sum(c*v for c, (_, v) in zip(coeffs, bvals))
        reds[(m, j)] = coeffs
        terms = ", ".join(f"{c}*{nm}" for c, (nm, _) in zip(coeffs, bvals) if c)
        print(f"I_{{{m},{j}}} (w{n}) = {terms}   |resid| {mp.nstr(abs(res),3)}", flush=True)
pickle.dump(reds, open('I_reductions_w2_7.pkl','wb'))
log("PSLQ reductions of I_{m,j} for weights 2..7 computed")
