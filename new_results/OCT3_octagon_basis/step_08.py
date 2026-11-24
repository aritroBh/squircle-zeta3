import json
from mpmath import mp, mpf, nstr
mp.dps = 280
d = json.load(open("consts62.json"))
X = [mpf(s) for s in d["vals"]]
# known relation: 24 Li3(1/2) - 21 z3 + 2 pi2 l2 - 4 l2^3 = 0 ; positions 22,34,36,40
s = 24*X[22] - 21*X[34] + 2*X[36] - 4*X[40]
print("direct residual of known relation:", nstr(s, 5))

B = 10**255
v = B*X[0]
print("type:", type(v))
r = round(v)
print("round type:", type(r), "-> int ok:", int(r))
print("check |X0 - round(BX0)/B|:", nstr(abs(X[0] - mpf(int(r))/B), 5))
