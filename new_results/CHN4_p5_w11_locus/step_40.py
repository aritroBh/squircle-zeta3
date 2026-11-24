import sympy as sp
# Peek at the structure of the files
for fn in ['R_a1a2a4.txt','R_factors.txt','subres1_A.txt','subres1_B.txt']:
    s=open(fn).read()
    print(f"{fn}: {len(s)} chars, head: {s[:120]!r}")
    print()
