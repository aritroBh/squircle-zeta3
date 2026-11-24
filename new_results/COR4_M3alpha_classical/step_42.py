import sys
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *

# simple test: shuffle of (e1) with a 3-word must have delta22 = 0
e1 = (Fr(0), Fr(1)); e0 = (Fr(1), Fr(0))
b3 = ((Fr(0),Fr(-1)),)*3
S = {}
for w in shuffles((e1,), b3):
    S[w] = (Fr(1), Fr(0))
print("shuffle words:", list(S))
d = delta22_spec(S, 2)
print("delta22(shuffle e1 x b^3):", d if d else "= 0 OK")

# and shuffle of 2-word with 2-word
S2 = {}
for w in shuffles((e0, e1), (e1, e0)):
    S2[w] = (Fr(1), Fr(0))
d2 = delta22_spec(S2, 2)
print("delta22(shuffle (01)x(10)):", d2 if d2 else "= 0 OK")

# delta22 of a single word (e1, b, b, b):
S3 = {(e1,) + b3: (Fr(1), Fr(0))}
print("delta22 single (e1,b,b,b):", delta22_spec(S3, 2))
