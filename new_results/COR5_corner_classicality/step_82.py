import sys, json
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import *

Q = W4Quotient(2)
# H_{0,1,0,1}(c): S = g(1)g(0)g(1)g(0) = (-[1-c])[c](-[1-c])[c]
# specialized: -[1-1/2] = +[2] = (1,0); [c] = -[2] = (-1,0)
word = ((Fr(1),Fr(0)), (Fr(-1),Fr(0)), (Fr(1),Fr(0)), (Fr(-1),Fr(0)))
print("quotient coords of S(H_{0,1,0,1}(1/2)):", Q.coords({word: Fr(1)}), "(expect (0,0,0): Li22(1/2,1) is products-only)")

# also H_{0,0,1,1}(1/2) = Li31(1/2,1) = -1/8 zeta3 log2 + pi4/720 + log2^4/24: also products-only
word2 = tuple([(Fr(1),Fr(0)), (Fr(1),Fr(0)), (Fr(-1),Fr(0)), (Fr(-1),Fr(0))])  # (-[1-c])(-[1-c])[c][c] spec
print("quotient coords of S(H_{0,0,1,1}(1/2)):", Q.coords({word2: Fr(1)}), "(expect (0,0,0): Li31(1/2,1) is products-only)")

# H_{0,0,0,1}(1/2) = Li4(1/2): projection (0,0,0) per li4v table
word3 = tuple([(Fr(1),Fr(0))] + [(Fr(-1),Fr(0))]*3)
print("quotient coords of S(H_{0,0,0,1}(1/2)):", Q.coords({word3: Fr(1)}), "(expect (0,0,0) per table)")
