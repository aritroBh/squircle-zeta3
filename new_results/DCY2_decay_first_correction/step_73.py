from mpmath import mp, gammainc
mp.dps = 40
for s,x in [(6.5, 1.9), (5.5, 1.9), (4.5, 1.9), (6.5, 0.742), (6.5, 1.8016)]:
    print(f"Gamma({s}, {x}) =", gammainc(s, x, regularized=False))
# and the difference relevant for the edge
print("diff:", gammainc(6.5, 0.742, regularized=False) - gammainc(6.5, 1.8016, regularized=False))
