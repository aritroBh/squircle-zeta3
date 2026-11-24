
from fractions import Fraction as Fr
C0 = Fr(1,2)
w = (1,0,1,0,1)
Delta1 = s05.delta_word_sym(w, 0, C0)
S5w = s05.S_J(w, C0, 0)
D5w = s05.delta_minus_5sym(S5w)
print("delta_word_sym:", len(Delta1), " delta^-(5sym):", len(D5w))
keys = sorted(set(Delta1) | set(D5w))
for k in keys:
    a = Delta1.get(k, Fr(0)); b = D5w.get(k, Fr(0))
    ratio = (a / b) if b else None
    print(f"  {k}: delta_word={a}  delta^-={b}  ratio={ratio}")
