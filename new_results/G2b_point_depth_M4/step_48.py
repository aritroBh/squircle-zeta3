
from fractions import Fraction as Fr
C0 = Fr(1,2)
Delta = {}
for w, coef in s05.signed_words():
    for lw, v in s05.delta_word_sym(w, 0, C0).items():
        Delta[lw] = Delta.get(lw, Fr(0)) + coef * v
Delta = {k: v for k, v in Delta.items() if v}
S5 = s05.omega5_half()
D5 = s05.delta_minus_5sym(S5)
print("Delta:", len(Delta), " D5:", len(D5))
keys = sorted(set(Delta) | set(D5))
ratios = set()
for k in keys:
    a = Delta.get(k, Fr(0)); b = D5.get(k, Fr(0))
    r = (a / b) if b else None
    ratios.add(r)
    if len(keys) < 40:
        print(f"  {k}: delta_word={a}  delta^-={b}  ratio={r}")
print("distinct ratios (delta_word / delta^-):", ratios)
