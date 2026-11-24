import sys
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
import cor4_sym as c4
import cor5_sym as cs

# square corner: delta_22(F(1/sqrt(2))) = 0 with letters in span{[2],[1+sqrt2]}
# specialize F_symbol.json: letter map old basis (2,c,1-c,1+c) -> new basis (2, 1+sqrt2):
# [c] = [1/sqrt2] = -(1/2)[2]; [1-c] = [1-1/sqrt2] = [(sqrt2-1)/sqrt2] = [b] - (1/2)[2] with b = sqrt2-1 = (1+sqrt2)^{-1}
# [1+c] = [(sqrt2+1)/sqrt2] = [1+sqrt2] - (1/2)[2]
base = "<project>/research"
Fterms = cs.parse_F_symbol(base + "/08_theory/COR4_M3alpha_classical/F_symbol.json")
# carry c^e = (1/sqrt2)^e = 2^{-e/2}: coefficients become half-integers in the exponent of 2 only:
# track scalar as (Fr exponent of 2^(-e/2)): since only the [2] letter carries it, fold into letters? No:
# delta22 ignores scalars except through coefficients: use tuple (Fr,Fr) = (rational, sqrt2-part)
H = Fr(1,2)
letter_map = [(Fr(1), Fr(0)), (Fr(-1,2), Fr(0)), (Fr(-1,2), Fr(-1)), (Fr(-1,2), Fr(1))]
# scalar c^e = 2^{-e/2}: represent as (Fr(-e,2) as power of 2 with half-integer exponent)
F_sq = {}
for e, word, cf in Fterms:
    nw = tuple(cs.lincomb(l, letter_map) for l in word)
    cval = Fr(2) ** (-e) if e % 2 == 0 else (Fr(2) ** (-(e+1)//2) if e < 0 else Fr(1, 2**( (e-1)//2 )) / 2)
    # c = 2^{-1/2}; c^e = 2^{-e/2}: e even -> 2^{-e/2} (Fr); e odd -> 2^{-e/2} = 2^{-(e+1)/2} * sqrt2 (sqrt2 factor)
    F_sq[nw] = F_sq.get(nw, Fr(0)) + cf * (Fr(2) ** (Fr(-e, 2)) if False else Fr(2) ** (-e) if e % 2 == 0 else Fr(0))
# that got messy; instead: scale each word by c^e exactly with (p,q) pairs meaning p + q*sqrt2
def cpair(e):
    # c^e with c = 2^{-1/2}: e even: (2^{-e//2}, 0); e odd: (2^{-(e+1)//2}, 0)*sqrt2 -> (0, 2^{-(e+1)//2})
    if e % 2 == 0:
        return (Fr(2) ** (-e // 2), Fr(0))
    else:
        return (Fr(0), Fr(2) ** (-(e + 1) // 2) if e < 0 else Fr(1, 2 ** ((e + 1) // 2)) * 2)
F_sq = {}
for e, word, cf in Fterms:
    nw = tuple(cs.lincomb(l, letter_map) for l in word)
    p, q = cpair(e)
    cur = F_sq.get(nw, (Fr(0), Fr(0)))
    F_sq[nw] = (cur[0] + cf * p, cur[1] + cf * q)
F_sq = {w: c for w, c in F_sq.items() if c != (Fr(0), Fr(0))}
print("S(F(1/sqrt2)): terms:", len(F_sq))
d22 = c4.delta22L2_spec(F_sq, 2)
print("delta_22(F(1/sqrt2)):", d22 if d22 else "0  (as required by thm:W4R)")

# thm:W4R's seven arguments are {sqrt2, sqrt2-1}-units:
sq2 = (Fr(0), Fr(1))
def cmul(a, b): return (a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
args = [("b", (Fr(-1), Fr(1))), ("-b", (Fr(1), Fr(-1))), ("-b^2", (Fr(3), Fr(-2))), ("1/sqrt2", (Fr(0), Fr(1,2))),
        ("1/2", (Fr(1,2), Fr(0))), ("1-1/sqrt2", (Fr(1), Fr(-1,2))), ("2-sqrt2", (Fr(2), Fr(-1)))]
for name, z in args:
    one_minus = (Fr(1) - z[0], Fr(0) - z[1])
    vx = c4.sunit_vec_sq2(z)
    vy = c4.sunit_vec_sq2(one_minus)
    print(f"arg {name}: z={z} sunit {vx}; 1-z={one_minus} sunit {vy}")
