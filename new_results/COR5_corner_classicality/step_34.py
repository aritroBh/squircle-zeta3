import sys, json
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *

base = "<project>/research"
cor4 = base + "/08_theory/COR4_M3alpha_classical"

# --- F's function symbol, specialized at c = 1/2 ---
# old basis: 0=[2], 1=[c], 2=[1-c], 3=[1+c]; new basis: 0=[2], 1=[3]
Fterms = parse_F_symbol(cor4 + "/F_symbol.json")
letter_map = [(Fr(1), Fr(0)),   # [2] -> [2]
              (Fr(-1), Fr(0)),  # [c]=[1/2] -> -[2]
              (Fr(-1), Fr(0)),  # [1-c]=[1/2] -> -[2]
              (Fr(-1), Fr(1))]  # [1+c]=[3/2] -> -[2]+[3]
F_half = {}
for e, word, cf in Fterms:
    cval = Fr(2) ** (-e)          # c^e at c=1/2 is 2^{-e}
    nw = tuple(lincomb(l, letter_map) for l in word)
    F_half[nw] = F_half.get(nw, Fr(0)) + cf * cval
F_half = {w: c for w, c in F_half.items() if c}
print("S(F(1/2)): terms:", len(F_half))
for w, c in sorted(F_half.items()):
    def ls(v):
        a, b = v
        s = ""
        if a: s += f"[2]^{a}" if a != 1 else "[2]"
        if b: s += ("*" if s else "") + (f"[3]^{b}" if b != 1 else "[3]")
        return s or "1"
    print(f"  {c!s:>5}  " + " ".join(f"{ls(v):<12}" for v in w))

Q = W4Quotient(2)
vF = Q.coords(F_half)
print()
print("pi(S(F(1/2))) in Q^4 coords:", vF)
print("quotient basis words:",
      ["".join('2' if l == Q.basis_letters[0] else '3' for l in Q.W4[j]) for j in Q.free])
