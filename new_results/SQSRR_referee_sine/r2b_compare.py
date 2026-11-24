import re, sys
from mpmath import mp, mpf, nstr
from law import law
mp.dps = 80
for fn in ('r2_defquad_A.log', 'r2_defquad_B.log'):
    txt = open(fn).read()
    for n, b1, b2 in re.findall(r'n (\d+) b\(2\) \S+\n\s+b1\s+(\S+)\n\s+b2\s+(\S+)', txt):
        n = int(n); be, lb2 = law(n)
        print(fn, 'n=%3d' % n, "b_n'' quad =", nstr(mpf(b2), 45), " quad-law: b'", nstr(mpf(b1)-be, 3), " b''", nstr(mpf(b2)-lb2, 3))
