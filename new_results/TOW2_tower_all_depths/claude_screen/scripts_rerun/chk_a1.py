import sympy as sp, tow2lib as T
d = T.d; i = sp.I
x = sp.radsimp(sp.simplify(T.complete_bell(1) + i*(d-2)))
print('L(1) + i(d-2) simplified:', x)
for r in range(1,5):
    print('L(%d) contains sqrt2 after radsimp:'%r, sp.radsimp(sp.simplify(T.L(r))).has(sp.sqrt(2)))
