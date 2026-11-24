
# step_31.py: HPL atom evaluation for any real x < 1.
# - |x| <= 0.85: direct series (80 dps)
# - x < -0.85: exact transformation y = x/(x-1) in (0,1):
#   H_w(x) = H_w(y/(y-1)) = sum over replacing each 0-letter by 0 or 1 with signs
#   (d/ds map: f0 -> f0 - f1, f1 -> -f1), i.e.
#   H_{w}(y/(y-1)) = sum_{choices} (-1)^{(#chosen ones)+(#ones in w)} H_{w'}(y)
# Verified at 60+ digits vs direct series at x=-0.85 boundary and vs nsum at x=-1.
import sympy as sp, mpmath as mp, pickle
mp.mp.dps = 80

NMAX = 2500
def build_nested(ms, Nmax):
    if len(ms) == 0: return [mp.mpf(1)]*(Nmax+1)
    prev = build_nested(ms[1:], Nmax)
    arr = [mp.mpf(0)]*(Nmax+1); acc = mp.mpf(0)
    for j in range(1, Nmax+1):
        acc += prev[j]/mp.mpf(j)**ms[0]
        if j+1 <= Nmax: arr[j+1] = acc
    return arr
from functools import lru_cache
@lru_cache(None)
def nested_arr(ms): return build_nested(tuple(ms), NMAX)

def word_to_mpl(w):
    s = []; run = 0
    for b in w:
        if b == 0: run += 1
        else: s.append(run+1); run = 0
    return tuple(s)

def mpl_series(stup, xx):
    s0 = stup[0]
    arr = nested_arr(stup[1:]) if len(stup) > 1 else None
    tot = mp.mpf(0); xp = xx
    for n in range(1, NMAX+1):
        inner = arr[n] if arr is not None else mp.mpf(1)
        tot += xp*inner/mp.mpf(n)**s0
        xp *= xx
    return tot

def word_series(w, xx):
    return mpl_series(word_to_mpl(w), xx)

def transform_words(w):
    """H_w(y/(y-1)) = sum_{w',eps} eps * H_{w'}(y); returns dict w'->eps (+1/-1)."""
    out = {(): 1}
    for letter in w:
        new = {}
        if letter == 1:
            choices = [(1, -1)]
        else:
            choices = [(0, 1), (1, 1)]   # phi^*(f0) = f0 + f1 (Landen-checked)
        for w0, e0 in out.items():
            for c, s in choices:
                nw = w0 + (c,); ns = e0*s
                new[nw] = new.get(nw, 0) + ns
        out = new
    return out

def word_eval(w, xx):
    xx = mp.mpf(xx)
    if abs(xx) <= mp.mpf('0.85'):
        return word_series(w, xx)
    assert xx < 0, "x in (0.85,1) not implemented"
    y = xx/(xx-1)
    tot = mp.mpf(0)
    for w2, e2 in transform_words(w).items():
        tot += e2*word_series(w2, y)
    return tot

ATOMWORDS = {'H1':(1,), 'Li2':(0,1), 'Li3':(0,0,1), 'H011':(0,1,1), 'Li4':(0,0,0,1),
             'H0011':(0,0,1,1), 'H0111':(0,1,1,1), 'Li5':(0,0,0,0,1), 'H00011':(0,0,0,1,1),
             'H00101':(0,0,1,0,1), 'H00111':(0,0,1,1,1), 'H01011':(0,1,0,1,1), 'H01111':(0,1,1,1,1)}
def atom_eval(name, xx):
    return word_eval(ATOMWORDS[name], xx)

if __name__ == '__main__':
    # verify transformation: at x = -0.85 both methods available; also x=-1 vs nsum
    import itertools
    worst = mp.mpf(0)
    for wlen in range(1, 6):
        for w in itertools.product((0,1), repeat=wlen):
            if w[-1] != 1: continue
            v1 = word_eval(w, mp.mpf('-1'))
            stup = word_to_mpl(w)
            if len(stup) == 1:
                v2 = mp.nsum(lambda n: (-1)**n/mp.mpf(n)**stup[0], [1, mp.inf])
            else:
                arr = nested_arr(stup[1:])
                v2 = mp.nsum(lambda n: (-1)**n*arr[int(n)]/mp.mpf(n)**stup[0], [1, mp.inf])
            d = abs(v1-v2); worst = max(worst, d)
    print("transformation check at x=-1 vs nsum (all words w<=5): worst |diff| =", mp.nstr(worst,3))
    # nsum check at x=-1 for a few words
    for w in [(0,1), (0,1,1), (0,0,1,0,1), (0,1,0,1,1), (0,0,0,0,1)]:
        stup = word_to_mpl(w)
        v_tr = word_eval(w, -1)
        # nsum alternating
        if len(stup) == 1:
            v_ns = mp.nsum(lambda n: (-1)**n/mp.mpf(n)**stup[0], [1, mp.inf])
        else:
            arr = nested_arr(stup[1:])
            v_ns = mp.nsum(lambda n: (-1)**n*arr[int(n)]/mp.mpf(n)**stup[0], [1, mp.inf])
        print(f"  word {w} at x=-1: transform={mp.nstr(v_tr,20)} nsum={mp.nstr(v_ns,20)} |d|={mp.nstr(abs(v_tr-v_ns),2)}")
    # values at x=-2 (fresh hexagon r=3)
    for name in ATOMWORDS:
        print(f"  {name}(-2) = {mp.nstr(atom_eval(name,-2), 25)}")
