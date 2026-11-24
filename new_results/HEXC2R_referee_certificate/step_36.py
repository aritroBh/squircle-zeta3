
reflib_numeric = '''

# ================= independent numeric evaluator =================
# Values from the DEFINITION as iterated sums. For a word with letters = forms
# (dt/t or a*dt/(1-a*t)), scanned into blocks ((s_1,a_1),...,(s_r,a_r)):
#     I(w; 0->x) = sum_{n_1>...>n_r>=1}  x^{n_1} prod_i a_i^{n_i-n_{i+1}} / prod_i n_i^{s_i}
# (n_{r+1}=0). With |a_i| <= 1 and |x| = 1/2 every summand is O(2^{-n_1} poly), so the
# sum is evaluated by truncation. Equivalently with z_1 = a_1 x, z_i = a_i/a_{i-1}:
#     value = sum_{n_2>=1} z_2^{n_2}/n_2^{s_2} * L1(n_2) * Q(n_2),
#     L1(n) = sum_{n_1>n} z_1^{n_1}/n_1^{s_1},  Q = inner multiple sum (blocks 3..r),
# all accumulated incrementally in O(K) operations, K ~ (dps+20)/log10(2).
import mpmath as mp

_STATE2 = {}

def _blocks(letters):
    """letters: tuple of ('T',) or ('A', alpha). Returns block tuple ((s,alpha),...)."""
    out = []; s = 1
    for L in letters:
        if L[0] == 'T':
            s += 1
        else:
            out.append((s, L[1])); s = 1
    if s != 1:
        raise ValueError("divergent word (trailing dt/t)")
    return tuple(out)

def _zpow_list(z, K):
    out = [mp.mpc(1)] * (K + 1)
    for i in range(1, K + 1): out[i] = out[i - 1] * z
    return out

def msum(blocks, x):
    """the truncated multiple sum, |a_i|<=1, |x|<=1/2."""
    r = len(blocks)
    if r == 0: return mp.mpc(1)
    dps = mp.mp.dps
    K = int((dps + 20) / 0.30103) + 30
    s = [b[0] for b in blocks]
    a = [mp.mpc(b[1]) for b in blocks]
    x = mp.mpc(x)
    z1 = a[0] * x
    z1p = _zpow_list(z1, K)
    L = [mp.mpc(0)] * (K + 2)
    for n1 in range(K, 0, -1):
        L[n1 - 1] = L[n1] + z1p[n1] / mp.mpf(n1) ** s[0]
    if r == 1:
        return L[0]
    z2 = a[1] / a[0]
    total = mp.mpc(0)
    if r == 2:
        z2p = _zpow_list(z2, K)
        for n2 in range(1, K):
            total += z2p[n2] / mp.mpf(n2) ** s[1] * L[n2]
        return total
    if r == 3:
        z3 = a[2] / a[1]
        z2p = _zpow_list(z2, K); z3p = _zpow_list(z3, K)
        Q = mp.mpc(0)
        for n2 in range(1, K):
            total += z2p[n2] / mp.mpf(n2) ** s[1] * L[n2] * Q
            Q += z3p[n2] / mp.mpf(n2) ** s[2]
        return total
    if r == 4:
        z3 = a[2] / a[1]; z4 = a[3] / a[2]
        z2p = _zpow_list(z2, K); z3p = _zpow_list(z3, K); z4p = _zpow_list(z4, K)
        P = mp.mpc(0); Q = mp.mpc(0)
        for n2 in range(1, K):
            total += z2p[n2] / mp.mpf(n2) ** s[1] * L[n2] * Q
            Q += z3p[n2] / mp.mpf(n2) ** s[2] * P
            P += z4p[n2] / mp.mpf(n2) ** s[3]
        return total
    raise ValueError("depth > 4 not needed here")

# ---- CMZV words (letters -1, 0..5) ----
def _zeta():
    return mp.e ** (mp.pi * mp.j / 3)

def _aforms():
    Z = _zeta()
    return [None] + [ -Z**j / (1 - Z**j) for j in range(1, 6) ]

def word_letters_A(w):
    """CMZV word -> form letters for the 0 -> x side."""
    Z = _zeta()
    return tuple(('T',) if j == -1 else ('A', Z**j) for j in w)

def word_letters_B(u):
    """CMZV word u (for the x -> 1 side) -> (sign, form letters) after s = 1-t reflection:
    I(u; 1/2->1) = (-1)^{|u|} I(ref(u); 0->1/2), ref = reversed word with f(t) -> f(1-t):
    letter -1: ds/(1-s) (A,1) ; letter 0: ds/s ('T') ; letter j: -a_j ds/(1-a_j s)."""
    AF = _aforms()
    sign = (-1) ** len(u)
    out = []
    for j in reversed(u):
        if j == -1:
            out.append(('A', mp.mpc(1)))
        elif j == 0:
            out.append(('T',))
        else:
            out.append(('A', AF[j])); sign = -sign
    return sign, tuple(out)

def cmzv_value(w):
    """regularized value of an ADMISSIBLE level-6 CMZV word (first letter outermost)."""
    w = tuple(w)
    if not w: return mp.mpc(1)
    assert w[0] != 0 and w[-1] != -1, "word not admissible"
    tot = mp.mpc(0)
    for k in range(len(w) + 1):
        pre, suf = w[:k], w[k:]
        A = msum(_blocks(word_letters_A(suf)), mp.mpf(1) / 2)
        if pre:
            sgn, letters = word_letters_B(pre)
            B = sgn * msum(_blocks(letters), mp.mpf(1) / 2)
        else:
            B = mp.mpc(1)
        tot += B * A
    return tot

def cmzv_combo(D):
    """D: dict word -> Fraction coefficient."""
    tot = mp.mpc(0)
    for w, c in D.items():
        tot += mp.mpf(c.numerator) / c.denominator * cmzv_value(w)
    return tot

# ---- HPL words (letters 0, +/-1 = dt/t, dt/(1-t), dt/(1+t)), first letter outermost ----
def hpl_value(w, x):
    w = tuple(w)
    if not w: return mp.mpc(1)
    letters = tuple(('T',) if j == 0 else ('A', mp.mpc(j)) for j in w)
    return msum(_blocks(letters), mp.mpf(x))
'''
with open("reflib.py", "a") as f:
    f.write(reflib_numeric)
print("appended numeric evaluator, total chars:", len(open("reflib.py").read()))
