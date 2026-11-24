# fix path6 in hexc2_common.py
src = open('hexc2_common.py').read()
old = """    out = {}
    for k in range(len(w) + 1):
        left = regularize({shift(tuple(reversed(w[:k])), c): Fr(1)})
        right = regularize({shift(w[k:], e): Fr(1)})
        out = add((out, 1), (mul(left, right), Fr((-1) ** k)))
    return regularize(out)"""
new = """    # I(w; a->e) = sum_k I(w[:k]; 0->e) * I(w[k:]; a->0),  I(v; a->0) = (-1)^|v| I(rev(v); 0->a)
    out = {}
    n = len(w)
    for k in range(n + 1):
        left = regularize({shift(w[:k], e): Fr(1)})                    # I(w[:k]; 0 -> zeta^e)
        right = regularize({shift(tuple(reversed(w[k:])), c): Fr(1)})  # I(rev(w[k:]); 0 -> zeta^c)
        out = add((out, 1), (mul(left, right), Fr((-1) ** (n - k))))
    return regularize(out)"""
assert old in src
open('hexc2_common.py','w').write(src.replace(old, new))
# fix cmzv_value in numeric_cmzv.py
src = open('numeric_cmzv.py').read()
old = "    return sum(A(w[:k]) * B(w[k:]) for k in range(len(w) + 1))"
new = "    return sum(B(w[:k]) * A(w[k:]) for k in range(len(w) + 1))"
assert old in src
open('numeric_cmzv.py','w').write(src.replace(old, new))
print('fixed both')
import sys; sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import numeric_cmzv as NC
import mpmath as mp
NC.set_dps(50)
tests = [((3,), -mp.log(2)), ((1,), mp.pi*mp.j/3), ((-1,0), mp.zeta(2)), ((-1,-1,0), mp.zeta(3)),
         ((-1,-1,-1,0), mp.zeta(4)), ((2,), -mp.log(3)/2 + mp.pi*mp.j/6),
         ((-1,0,0), mp.polylog(2,1)*0 + mp.mpf('0.5822405264650125')),  # placeholder
        ]
for w, expect in tests[:6]:
    got = NC.cmzv_value(w); err = abs(got-expect)
    print(w, 'err=%.2e'%float(err), 'OK' if err < mp.mpf(10)**(-45) else 'FAIL')
# depth-2 coloured checks vs mpmath polylog: word (-1,-1,j) = Li_3(zeta^j)
for j in range(6):
    got = NC.cmzv_value((-1,-1,j)); expect = mp.polylog(3, mp.e**(mp.pi*mp.j*j/3))
    err = abs(got-expect)
    print((-1,-1,j), 'err=%.2e'%float(err), 'OK' if err < mp.mpf(10)**(-45) else 'FAIL', got if err>1e-40 else '')
