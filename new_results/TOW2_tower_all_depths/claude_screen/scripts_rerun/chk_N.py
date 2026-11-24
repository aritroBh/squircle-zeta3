import re, sympy as sp, tow2lib as T
d, m = T.d, T.m; i = sp.I
ns = {'sp': sp, 'd': d, 'm': m, 'i': i}
def getN(fn):
    txt = open(fn).read()
    j = txt.index('N = {')
    # find matching brace
    depth = 0; k = j + 4
    while True:
        c = txt[k]
        if c == '{': depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0: break
        k += 1
    return eval(txt[j+4:k+1], ns)
N12 = getN('prove_S12.py'); N34 = getN('prove_S34.py')
for l, Nl in {**N12, **N34}.items():
    print('N', l, 'equal to tow2lib.N:', sp.expand(Nl - T.N[l]) == 0)
