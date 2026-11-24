"""A=B=1: compare n*(f_n/law - 1) with the paper's (3 - 3 pi (-1)^n)/8 (main.tex, after thm:zero)."""
import os, pickle
from mpmath import mp, mpf, sqrt, cos, pi, factorial, nstr
mp.dps = 50
I = [mpf(x) for x in pickle.load(open(os.path.join(os.path.dirname(__file__), "In_direct.pkl"), "rb"))[(1, 1)]]
for n in range(30, 41):
    law = sqrt(2) * cos(pi * (3 - 2 * n) / 4) * factorial(n) * (2 / pi) ** n * mpf(n) ** -1.5
    print(n, nstr(n * (I[n] / law - 1), 6), nstr((3 - 3 * pi * (-1) ** n) / 8, 6))
