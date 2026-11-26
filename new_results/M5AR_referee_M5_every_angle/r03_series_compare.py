"""r03 (referee, own code): exact comparison of the Taylor/Laurent coefficients of every claimed part
(term lists c^e He[w] := c^e H_w(c) + (-c)^e H_w(-c), H_w from its definition) with the independent exact
mu_{k,n} of r01 (Beta/polygamma route), for n <= NMAX; also the closed forms of result.md and of the paper.
Does not import m5fit, g8_allN, g6_common."""
import json, re, sys, pickle
from fractions import Fraction as Fr
from functools import lru_cache

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
MU = pickle.load(open(sys.argv[2] if len(sys.argv) > 2 else "r01_mu_k0to5_n150.pkl", "rb"))
NT = 2*NMAX + 8          # highest power of c computed
LO = -6

@lru_cache(maxsize=None)
def hcoef(w):
    """tuple a[0..NT]: [x^m] H_w(x), H_1=-log(1-x), H_-1=log(1+x), H_{a,w}=int_0^x f_a H_w, f_0=1/t, f_1=1/(1-t), f_-1=1/(1+t)."""
    if len(w) == 1:
        a = w[0]; assert a in (1, -1)
        return tuple([Fr(0)] + [Fr(1, m) if a == 1 else Fr((-1)**(m+1), m) for m in range(1, NT+1)])
    inner = hcoef(w[1:]); a = w[0]
    out = [Fr(0)]*(NT+1)
    if a == 0:
        for m in range(1, NT+1): out[m] = inner[m]/m
    else:
        # f_a(t) H_w(t) = sum_j inner[j] t^j * sum_l a^l t^l  (f_1 = 1/(1-t): a=1 ; f_-1 = 1/(1+t): a=-1)
        s = Fr(0)
        for m in range(1, NT+1):
            s = s*a + inner[m-1]       # [t^(m-1)] of f_a H_w
            out[m] = s/m
    return tuple(out)

def series_of_terms(raw):
    """Laurent coefficients {N: Fraction}, LO <= N <= NT-? , of the term list."""
    C = {N: Fr(0) for N in range(LO, NT - 4)}
    for nm, v in raw:
        q = Fr(v)
        if nm == '1': C[0] += q; continue
        m = re.match(r'c\^(-?\d+)$', nm)
        if m: C[int(m.group(1))] += q; continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]$', nm)
        e = int(m.group(1)); w = tuple(int(t) for t in m.group(3).split(','))
        a = hcoef(w)
        for N in C:
            mm = N - e
            if 1 <= mm <= NT:
                x = a[mm]
                if m.group(2) == 'He':
                    C[N] += q*(x + (x if (e + mm) % 2 == 0 else -x))      # c^e H_w(c) + (-c)^e H_w(-c)  (integer parity; (-1)**e is a float for e<0)
                else:
                    C[N] += q*x
    return C

def from_func(f):
    """series {N: Fr} from a dict {N: coefficient}."""
    return {N: f.get(N, Fr(0)) for N in range(LO, NT - 4)}

def lin(*pairs):
    out = {}
    for a, S in pairs:
        for N, v in S.items(): out[N] = out.get(N, 0) + a*v
    return out

def mu_target(k, mon):
    """series sum_n [mon] mu_{k,n} c^{2n}."""
    S = {N: Fr(0) for N in range(LO, NT - 4)}
    for n in range(1, NMAX+1):
        S[2*n] = MU[k][n].get(mon, Fr(0))
    return S

def compare(S, T, label):
    bad = [N for N in range(LO, 2*NMAX + 2) if S.get(N, 0) != T.get(N, 0)]
    print(f"  {label}: coefficients c^N, {LO} <= N <= {2*NMAX+1}: {'ALL EQUAL' if not bad else 'MISMATCH at N=' + str(bad[:8])}")
    return not bad

Z = lambda **kw: tuple(kw.get(f"z{w}", 0) for w in range(2, 6))
load = lambda f: json.load(open(f))
D = "claim_copy/"
ok = []
print(f"=== (A) each fitted part of M_5 vs independent mu_(5,n), n <= {NMAX}")
parts = {'1': Z(), 'z2': Z(z2=1), 'z2^2': Z(z2=2), 'z3': Z(z3=1), 'z4': Z(z4=1), 'z5': Z(z5=1), 'z2*z3': Z(z2=1, z3=1)}
S5 = {}
for nm, mon in parts.items():
    raw = load(D + f"hpl_fit_k5_{nm}.json")
    S5[nm] = series_of_terms(raw)
    ok.append(compare(S5[nm], mu_target(5, mon), f"[{nm}]M5 ({len(raw)} lines)"))
allmons = set(m for n in MU[5] for m in MU[5][n] if n <= NMAX)
print("  zeta monomials occurring in mu_(5,n), n<=NMAX:", sorted(allmons), " -> exactly the 7 parts:", allmons == set(parts.values()))

print("=== (B) closed forms [z5]mu5 = -3/((2n-1)(2n+1)(2n+3)), [z2z3]mu5 = 2n/((2n-1)(2n+1)(2n+3))")
ok.append(all(MU[5][n].get(Z(z5=1)) == Fr(-3, (2*n-1)*(2*n+1)*(2*n+3)) for n in range(1, NMAX+1)))
ok.append(all(MU[5][n].get(Z(z2=1, z3=1)) == Fr(2*n, (2*n-1)*(2*n+1)*(2*n+3)) for n in range(1, NMAX+1)))
print("  z5:", ok[-2], "  z2z3:", ok[-1])
P4 = lambda mon: mu_target(4, (mon + (0,))[:4])
ok.append(all(MU[4][n].get(Z(z2=2), 0) == Fr(-n, (2*n-1)*(2*n+1)*(2*n+3)) for n in range(1, NMAX+1)))
print("  [z2^2]mu_(4,n) = -n/((2n-1)(2n+1)(2n+3)) (result.md s6):", ok[-1])

print("=== (C) the ladder parts of the THEOREM as stated, with the paper's Phi, G8's Psi, Psi_2, and F_2")
# Phi = sum c^2n/(4n^2-1)
Phi = from_func({2*n: Fr(1, 4*n*n-1) for n in range(1, NT)})
# artanh c = sum_{m odd} c^m/m
def atanh_times(coeffs):   # coeffs {e: q} -> series of (sum q c^e) * artanh c
    out = {}
    for e, q in coeffs.items():
        for m in range(1, NT+1, 2):
            out[m+e] = out.get(m+e, 0) + q*Fr(1, m)
    return out
F2 = lin((1, from_func(atanh_times({1: Fr(-1, 16), -1: Fr(-1, 8), -3: Fr(3, 16)}))), (1, from_func({0: Fr(1, 16), -2: Fr(-3, 16)})))
ok.append(compare(F2, mu_target(4, Z(z2=2)), "F_2 := (-c/16-1/(8c)+3/(16c^3))artanh c+1/16-3/(16c^2)  vs [z2^2]M4 (own mu)"))
ok.append(compare(Phi, mu_target(2, Z(z2=1)), "Phi vs [z2]M2 (own mu)"))
Psi = series_of_terms(load(D + "g8_terms_k3_partz2.json"))
Psi2 = series_of_terms(load(D + "g8_terms_k4_partz2.json"))
ok.append(compare(Psi, mu_target(3, Z(z2=1)), "Psi (G8 terms) vs [z2]M3 (own mu)"))
ok.append(compare(Psi2, mu_target(4, Z(z2=1)), "Psi_2 (G8 terms) vs [z2]M4 (own mu)"))
G = S5['z2^2']
ok.append(compare(lin((-1, Psi2)), mu_target(5, Z(z3=1)), "-Psi_2 vs [z3]M5"))
ok.append(compare(lin((1, Psi), (1, G)), mu_target(5, Z(z4=1)), "Psi + G vs [z4]M5"))
ok.append(compare(lin((-1, Phi), (-2, F2)), mu_target(5, Z(z5=1)), "-Phi - 2F_2 vs [z5]M5"))
ok.append(compare(lin((-2, F2)), mu_target(5, Z(z2=1, z3=1)), "-2F_2 vs [z2z3]M5"))
# paper example rung: [z5]M_k = -[z2]M_{k-3} - 2[z2^2]M_{k-1} at k=5, and [z3]M_{k+1} = -[z2]M_k for k<=4
ok.append(all(MU[5][n].get(Z(z5=1), 0) == -MU[2][n].get(Z(z2=1), 0) - 2*MU[4][n].get(Z(z2=2), 0) for n in range(1, NMAX+1)))
print("  paper's example rung [z5]M5 = -[z2]M2 - 2[z2^2]M4 on own mu:", ok[-1])
ok.append(all(MU[k+1][n].get(Z(z3=1), 0) == -MU[k][n].get(Z(z2=1), 0) for k in range(2, 5) for n in range(1, NMAX+1)))
print("  [z3]M_{k+1} = -[z2]M_k, k=2,3,4, on own mu:", ok[-1])

print("=== (D) closed forms in result.md")
L = {2*m: Fr(-1, m) for m in range(1, NT)}                                       # log(1-c^2), untruncated before the c^-4 shift
def mul_laurent(coeffs, S):
    out = {}
    for e, q in coeffs.items():
        for N, v in (S.items() if isinstance(S, dict) else S):
            if v: out[N+e] = out.get(N+e, 0) + q*v
    return from_func(out)
A = from_func(atanh_times({0: Fr(1)}))
def Hser(w):
    a = hcoef(w); return from_func({m: a[m] for m in range(1, NT+1)})
Dd = lin((1, Hser((0, 1))), (1, Hser((0, -1))), (-1, Hser((-1, 1))), (1, Hser((1, -1))))
# G = 7/16 - 3/(4c^2) - (16-19c^2+c^4)/(32c^4) L + (7-5c^2)/(16c^3) A + (c^2+3)(c^2-1)/(32c^3) D
Gcf = lin((1, from_func({0: Fr(7, 16), -2: Fr(-3, 4)})),
          (1, mul_laurent({-4: Fr(-16, 32), -2: Fr(19, 32), 0: Fr(-1, 32)}, L)),
          (1, mul_laurent({-3: Fr(7, 16), -1: Fr(-5, 16)}, A)),
          (1, mul_laurent({1: Fr(1, 32), -1: Fr(2, 32), -3: Fr(-3, 32)}, Dd)))
ok.append(compare(Gcf, mu_target(5, Z(z2=2)), "G closed form vs [z2^2]M5 (own mu)"))
z5cf = lin((1, from_func({0: Fr(-5, 8), -2: Fr(3, 8)})), (1, mul_laurent({-3: Fr(-3, 8), -1: Fr(6, 8), 1: Fr(-3, 8)}, A)))
ok.append(compare(z5cf, mu_target(5, Z(z5=1)), "[z5] closed form -5/8+3/(8c^2)-3(1-c^2)^2/(8c^3)A"))
z23cf = lin((1, from_func({0: Fr(-1, 8), -2: Fr(3, 8)})), (1, mul_laurent({-3: Fr(-3, 8), -1: Fr(2, 8), 1: Fr(1, 8)}, A)))
ok.append(compare(z23cf, mu_target(5, Z(z2=1, z3=1)), "[z2z3] closed form -1/8+3/(8c^2)-(1-c^2)(3+c^2)/(8c^3)A"))
Psicf = lin((1, from_func({0: Fr(1, 2)})), (1, mul_laurent({-2: Fr(1, 2), 0: Fr(-3, 4)}, L)), (1, mul_laurent({-1: Fr(-1, 2)}, A)),
            (1, mul_laurent({-1: Fr(1, 4), 1: Fr(-1, 4)}, Dd)))
ok.append(compare(Psicf, mu_target(3, Z(z2=1)), "Psi closed form vs [z2]M3 (own mu)"))
# pi^4 part: zeta4 [z4] + zeta2^2 [z2^2] = pi^4 ([z4]/90 + [z2^2]/36)
pi4 = lin((Fr(1, 90), mu_target(5, Z(z4=1))), (Fr(1, 36), mu_target(5, Z(z2=2))))
pi4cf = lin((1, from_func({0: Fr(65, 2880), -2: Fr(-84, 2880)})),
            (1, mul_laurent({0: Fr(-55, 5760), -2: Fr(165, 5760), -4: Fr(-112, 5760)}, L)),
            (1, mul_laurent({-1: Fr(-51, 2880), -3: Fr(49, 2880)}, A)),
            (1, mul_laurent({1: Fr(-3, 1920), -1: Fr(10, 1920), -3: Fr(-7, 1920)}, Dd)))
ok.append(compare(pi4cf, pi4, "pi^4 bracket of result.md vs [z4]/90+[z2^2]/36 (own mu)"))
# paper Thm M4alpha pi^4 part: pi^4/2880 [(9c - 30/c + 21/c^3) artanh c + 23 - 21/c^2]
pi4_4 = lin((Fr(1, 90), mu_target(4, Z(z4=1))), (Fr(1, 36), mu_target(4, Z(z2=2))))
pi4_4cf = lin((1, mul_laurent({1: Fr(9, 2880), -1: Fr(-30, 2880), -3: Fr(21, 2880)}, A)), (1, from_func({0: Fr(23, 2880), -2: Fr(-21, 2880)})))
ok.append(compare(pi4_4cf, pi4_4, "paper Thm M4alpha pi^4 bracket vs own mu_4"))
print("D coefficient of the pi^4 bracket: (1-c^2)(3c^2-7)/1920 /c^3 = (-3c^4+10c^2-7)/(1920c^3): nonzero polynomial")
print("ALL CHECKS:", all(ok), f"({len(ok)} checks)")
