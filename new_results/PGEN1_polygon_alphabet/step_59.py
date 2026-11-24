# step_33: verify symbol recursion (B): S(I(0;w;1)) = sum_i S(I(remove a_i)) x v((a_{i+1}-a_i)/(a_i-a_{i-1}))
# chain (a_0=0, a_1..a_n, a_{n+1}=1). Numeric checks of G values vs symbol-consistent closed forms.
import mpmath as mp
mp.mp.dps = 50
def G_int(word, z=1, dps=50):
    """iterated integral G(a1..an; z) = int_0^z G(a1..a_{n-1};t) dt/(t-a_n) -- convention: outer letter last? 
    We use I(0;a1..an;1) = int_0^1 I(0;a1..a_{n-1};t) dt/(t-an)."""
    mp.mp.dps = dps
    if len(word) == 0: return mp.mpf(1)
    a = word[-1]
    sub = word[:-1]
    pts = [mp.mpf(0)]
    # split at singularities inside (0,1): poles where t = a in (0,1)
    for s in [a]:
        if 0 < s < 1: pts.append(mp.mpf(s))
    pts.append(mp.mpf(z))
    pts = sorted(set(pts))
    f = lambda t: G_int(sub, t, dps)/(t-a)
    return mp.quad(f, pts)

# test words: value of I(0;2;1) = -log2; I(0;3;1) = log(2/3)
print(G_int((2,)), -mp.log(2))
print(G_int((3,)), mp.log(mp.mpf(2)/3))
# G(2,3;1): predicted symbol -2 v2^v2 + v3^v2  (recursion B)
g23 = G_int((2,3))
cand = -mp.log(2)**2 + mp.log(2)*mp.log(3) + mp.re(mp.polylog(2,3))  # S(Li2(3)) = -v2(x)v3
print("G(2,3;1) =", g23)
print("candidate  =", cand)
print("diff:", g23-cand)
