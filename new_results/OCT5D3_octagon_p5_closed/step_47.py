import pickle, mpmath as mp
mp.mp.dps = 60
av = pickle.load(open("<project>/research/08_theory/OCT5D2_octagon_depth_test/atom_vals430.pkl","rb"))
print("OCT5D2 L32_b1:", av['L32_b1'][:60])
print("OCT5D2 L32_h1:", av['L32_h1'][:60])
# their function from step_76
def Li32_d2(x, y, N):
    tot = mp.mpf(0); xp = x
    for n in range(2, N+1):
        if n > 2: xp *= x
        inner = mp.mpf(0); yp = y
        for m in range(1, n):
            if m > 1: yp *= y
            inner += yp/(m*m)
        tot += xp/(n**3)*inner
    return tot
r2 = mp.sqrt(2); b = r2-1
print("their func (b,1,1600):", mp.nstr(Li32_d2(b, mp.mpf(1), 1600), 50))
print("mine                  :", mp.nstr(mp.mpf('0.02554867181188595195138362416354460326083'), 50))
# maybe their atom used x=1/b? or Li32(1, x)? or sum_{n>=m}?
def Li32_alt(x, y, N):  # n>=m variant
    tot = mp.mpf(0)
    for n in range(1, N+1):
        inner = mp.mpf(0); yp = mp.mpf(1)
        for m in range(1, n+1):
            if m > 1: yp *= y
            inner += yp/(m*m)
        tot += x**n/(n**3)*inner
    return tot
print("alt n>=m (b,1):", mp.nstr(Li32_alt(b, mp.mpf(1), 800), 50))
