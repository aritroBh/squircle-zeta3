
import sys, time, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import rul1, thx1_cert
import sympy as sp
import math

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

e_rm = pickle.load(open('e_rm_ext.pkl','rb'))
def E(r, m): return Fr(e_rm.get(f"{r},{m}", '0'))
def I_word_combo(m, j):
    return {(0,)*(m+1) + (-1,)*j: Fr(math.factorial(j)*math.factorial(m))}
def Dtilde_combo(n, m):
    out = {}
    for j in range(1, n-2*m):
        e = E(m, n-1-j)
        if e == 0: continue
        c = 2*e/math.factorial(j)
        for w, vv in I_word_combo(n-1-j, j).items():
            out[w] = out.get(w, Fr(0)) + c*vv
    return {k: v for k, v in out.items() if v}
def zeta_word(n): return {(0,)*(n-1)+(1,): Fr(-1)}
def zprod(*ns):
    out = {(): Fr(1)}
    for n in ns: out = rul1.sh_mul(out, zeta_word(n))
    return {k: v for k, v in out.items() if v}
Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))
def zbasis_combo(n):
    if n == 2: return [zprod(2)]
    if n == 3: return [zprod(3)]
    if n == 4: return [zprod(4)]
    if n == 5: return [zprod(5), zprod(2,3)]
    if n == 6: return [zprod(6), zprod(3,3)]
def V_target(n, m):
    out = Dtilde_combo(n, m)
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    for c, zb in zip(coeffs, zbasis_combo(n)):
        for w, vv in zb.items():
            out[w] = out.get(w, Fr(0)) - c*vv
    return {k: v for k, v in out.items() if v}

primes = [int(sp.prevprime(2**62)), int(sp.prevprime(2**61)), int(sp.prevprime(2**60))]
for W in [5, 6]:
    rels_data = rul1.relations_at_weight(W, families=('DS','DUAL'))
    rels = [r for fam, tag, r in rels_data]
    cols = sorted(rul1.conv_words_of_weight(W))
    st = pickle.load(open(f'struct_w{W}.pkl','rb'))
    for m in range(W//2):
        V = V_target(W, m)
        h, ok = thx1_cert.certify(f'D{W},{m}', V, rels, st, cols, primes)
        print(f"D~_{{{W},{m}}}: {'PASS' if ok else 'FAIL'}", flush=True)
        if ok:
            pickle.dump(h, open(f'cert_D{W}_{m}.pkl','wb'))
            log("D~_{%d,%d} certificate saved" % (W, m))

# ---- Q-polynomial assembly check: Q_n(u) = sum_m D~_{n,m} (-4u)^m/(2m)!  vs Qpolys.pkl ----
Qpolys = pickle.load(open('Qpolys.pkl','rb'))
u = sp.symbols('u')
z = {i: sp.Symbol(f'z{i}') for i in range(2, 10)}
zmap = {}
for n in [2,3,4,5,6,7,8,9]:
    if n in (2,3,4): zmap[n] = [z[n]]
    if n == 5: zmap[n] = [z[5], z[2]*z[3]]
    if n == 6: zmap[n] = [z[6], z[3]**2]
    if n == 7: zmap[n] = [z[7], z[2]*z[5], z[2]**2*z[3]]
    if n == 8: zmap[n] = [z[8], z[3]*z[5], z[2]*z[3]**2]
    if n == 9: zmap[n] = [z[9], z[2]*z[7], z[2]**2*z[5], z[2]**3*z[3], z[3]**3]
allok = True
for n in range(2, 10):
    Q = 0
    for m in range(n//2):
        coeffs = [sp.Rational(c) for c in Dtarget[f"{n},{m}"]]
        Dnm = sum(c*b for c, b in zip(coeffs, zmap[n]))
        Q += Dnm * (-4*u)**m / sp.factorial(2*m)
    Q = sp.expand(Q)
    target = sp.sympify(Qpolys[n], locals={f'z{i}': z[i] for i in range(2,10)} | {'u': u})
    diff = sp.expand(Q - target)
    ok = (diff == 0)
    allok &= ok
    print(f"Q_{n} assembly vs Qpolys.pkl: {'MATCH' if ok else 'MISMATCH: '+str(diff)}")
log("Q_n(u) assembly from D~ values matches Qpolys.pkl for all n<=9: %s" % allok)
