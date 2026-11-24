code = '''

# ================= target vectors (exact) =================
from fractions import Fraction as Fr

def prod(*vecs):
    out = {(): Fr(1)}
    for v in vecs: out = mul(out, v)
    return out

def Li(n, z): return li_atom(n, Fr(z))

def F4_hpl_vec():
    """weight-4 part of R3(1/2): sum cf * (1/2)^e * H_w(1/2)."""
    import json, re
    G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
    D = json.load(open(G1 + "/hpl_fit_k3_(0,0,0).json"))
    out = {}
    for key, val in D:
        m = re.match(r"c\\^(-?\\d+) H\\[(.*)\\]", key)
        if not m: continue
        e = int(m.group(1))
        w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        if len(w) != 4: continue
        out = add((out, 1), (hpl_to_cmzv(w), Fr(val) * Fr(1, 2) ** e))
    return out

def fit_vec(name, weight):
    """sum cf * (1/2)^e * H_w(1/2) over terms of given HPL weight in a G1 fit file."""
    import json, re
    G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
    D = json.load(open(G1 + "/" + name))
    out = {}
    for key, val in D:
        m = re.match(r"c\\^(-?\\d+) H\\[(.*)\\]", key)
        if key == "1": e, w = 0, ()
        elif key == "c": e, w = 1, ()
        else:
            e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        if len(w) != weight: continue
        v = {(): Fr(1)} if not w else hpl_to_cmzv(w)
        out = add((out, 1), (v, Fr(val) * Fr(1, 2) ** e))
    return out

def F4_closed_vec():
    L2, L3 = atom_log2(), atom_log3()
    Z3, PI2, PI4 = atom_zeta3(), smul(Fr(-1), atom_pi2()), atom_pi4()
    Li21, Li31, Li41 = Li(2, '1/3'), Li(3, '1/3'), Li(4, '1/3')
    Li3m, Li4m = Li(3, '-1/2'), Li(4, '-1/2')
    out = {}
    def T(vec, coef): 
        nonlocal out; out = add((out, 1), (vec, Fr(coef)))
    T(Li(4, -2), -21); T(Li41, Fr(-65,8)); T(Li(4, '2/3'), Fr(13,2)); T(Li(4, '1/2'), Fr(137,6))
    T(Li41, Fr(39,8)); T(Li(4, '-1/3'), Fr(39,8))         # 39/64 Li4(1/9) = 39/8 [Li4(1/3)+Li4(-1/3)]
    T(prod(Li31, L2), Fr(15,2)); T(prod(Li3m, L2), Fr(37,2)); T(prod(Li3m, L3), Fr(-5,4))
    T(prod(Z3, L2), Fr(119,16)); T(prod(Z3, L3), Fr(-37,32))
    T(prod(Li21, Li21), 1); T(prod(Li21, L2, L2), Fr(-23,4)); T(prod(Li21, L2, L3), Fr(7,4))
    T(prod(Li21, L3, L3), 1); T(prod(PI2, Li21), Fr(-1,6))
    T(prod(PI2, L2, L2), Fr(-455,144)); T(prod(PI2, L2, L3), Fr(9,8)); T(prod(PI2, L3, L3), Fr(-7,32))
    T(PI4, Fr(-949,2160))
    T(prod(L2,L2,L2,L2), Fr(-211,144)); T(prod(L2,L2,L2,L3), Fr(101,24)); T(prod(L2,L2,L3,L3), Fr(-15,8))
    T(prod(L2,L3,L3,L3), Fr(-35,24)); T(prod(L3,L3,L3,L3), Fr(113,192))
    return out

def W3_closed_vec():
    L2, L3 = atom_log2(), atom_log3()
    Z3, PI2 = atom_zeta3(), smul(Fr(-1), atom_pi2())
    Li21, Li31, Li3m = Li(2, '1/3'), Li(3, '1/3'), Li(3, '-1/2')
    out = {}
    def T(vec, coef):
        nonlocal out; out = add((out, 1), (vec, Fr(coef)))
    T(Z3, Fr(-15,4)); T(Li31, -3); T(Li3m, -8)
    T(prod(PI2, L2), Fr(-1,12)); T(prod(Li21, L2), 4); T(prod(Li21, L3), Fr(-3,2))
    T(prod(L2,L2,L2), Fr(1,3)); T(prod(L2,L2,L3), Fr(-5,2)); T(prod(L2,L3,L3), 2); T(prod(L3,L3,L3), Fr(-1,4))
    return out

def W2_closed_vec():
    L2, L32 = atom_log2(), atom_log32()
    return vadd(prod(L2, L2), Li(2, '1/3'), prod(L32, L32)) if False else \
        vadd(smul(Fr(1,2), prod(L2, L2)), Li(2, '1/3'), smul(Fr(1,2), prod(L32, L32)))

def target4():
    return vadd(F4_hpl_vec(), smul(Fr(-1), F4_closed_vec()))

def target3():
    # R3_3 + zeta(2) Psi1 - zeta(3)/2 - W3_closed
    R33 = fit_vec("hpl_fit_k3_(0,0,0).json", 3)
    Psi1 = fit_vec("hpl_fit_k3_(1,0,0).json", 1)
    lhs = vadd(R33, mul(atom_zeta2(), Psi1), smul(Fr(-1,2), atom_zeta3()))
    return vadd(lhs, smul(Fr(-1), W3_closed_vec()))

def target2():
    # R3_2 + zeta(2)/2 - W2_closed
    R32 = fit_vec("hpl_fit_k3_(0,0,0).json", 2)
    lhs = vadd(R32, smul(Fr(1,2), atom_zeta2()))
    return vadd(lhs, smul(Fr(-1), W2_closed_vec()))
'''
src = open('hexc2_common.py').read()
assert 'def target4' not in src
open('hexc2_common.py','a').write(code)
print('appended')
import sys; sys.dont_write_bytecode=True; sys.path.insert(0,'.')
from hexc2_common import target2, target3, target4
import mpmath as mp
import numeric_cmzv as NC
NC.set_dps(70)
def cmzv_combo(d):
    tot = mp.mpc(0)
    for w,c in d.items():
        tot += mp.mpf(c.numerator)/mp.mpf(c.denominator) * NC.cmzv_value(w)
    return tot
for name, T in [('T2', target2), ('T3', target3), ('T4', target4)]:
    v = T()
    val = cmzv_combo(v)
    print(name, 'support', len(v), ' value =', mp.nstr(val, 8))
