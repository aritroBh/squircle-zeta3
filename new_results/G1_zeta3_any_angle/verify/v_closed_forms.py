"""v_closed_forms.py -- verifier's re-typing (by hand, from result.md) of the claimed closed forms, evaluated
with mpmath polylog / own HPL power series, against the verifier's independent Cauchy-contour values
(runB/vc_*.txt, accurate to ~1e-46) at 13 angles (9 solver angles + 4 fresh: 0.3, 0.9, 2.6, 3.0 rad).

Checked: M0 (B0), M1 (B1), Phi closed form, W2 = (1/2)z2 - (1/2)Li2(c^2) (Thm C) both as the t-integral and as
the closed form, W3 x-integral table (Thm C), B2' (uniform classical), B2 (HPL form, own HPL series),
special-angle formula at 2pi/3, and the pi/3 Li3-table.
"""
import mpmath as mp

mp.mp.dps = 60
z2, z3 = mp.zeta(2), mp.zeta(3)
lam = mp.log(2)
ANG = {
    'pi_2': mp.pi/2, '2pi_3': 2*mp.pi/3, 'pi_3': mp.pi/3, '3pi_4': 3*mp.pi/4, 'pi_6': mp.pi/6,
    '5pi_6': 5*mp.pi/6, 'pi_4': mp.pi/4, 'gen1p3': mp.mpf(13)/10, 'gen2p0': mp.mpf(2),
    'new0p3': mp.mpf(3)/10, 'new0p9': mp.mpf(9)/10, 'new2p6': mp.mpf(26)/10, 'new3p0': mp.mpf(3),
}


ANG_hi = {'pi_2': lambda: mp.pi/2, '3pi_4': lambda: 3*mp.pi/4, 'pi_6': lambda: mp.pi/6, '5pi_6': lambda: 5*mp.pi/6,
          'pi_4': lambda: mp.pi/4, 'gen1p3': lambda: mp.mpf(13)/10, 'gen2p0': lambda: mp.mpf(2)}


def rd(fn):
    d = {}
    for l in open(fn):
        if l.startswith('M'):
            k, v = l.split('=')
            d[int(k.strip()[1])] = mp.mpf(v)
    return d


Li2 = lambda x: mp.polylog(2, x)
Li3 = lambda x: mp.polylog(3, x)


def M0f(c):
    return -(1 + c)*mp.log(1 + c) - (1 - c)*mp.log(1 - c)


def M1f(c):
    s = 0
    for e in (1, -1):
        s += (1 - e*c)*(Li2(e*c) + Li2((1 - e*c)/2))
    return z2 - lam**2 - lam*M0f(c) - s


def Phi_series(c):
    return mp.nsum(lambda n: c**(2*n)/(4*n**2 - 1), [1, mp.inf])


def Phi(c):
    return mp.mpf(1)/2 - (1 - c**2)/(2*c)*mp.atanh(c)


def M2_classical(c):
    lp, lm, l0 = mp.log(1 + c), mp.log(1 - c), mp.log(c)
    q = (1 - c**2)/(2*c)
    r = (1 - c)/(1 + c)
    v = z2*Phi(c) - Li2(c) - Li2(-c)
    v += (1 + c)*Li3(c) + (1 - c)*Li3(-c) + (1 - c)/c*Li3(1 - c) - (1 + c)/c*Li3(1/(1 + c)) - 2*c*Li3((1 + c)/2)
    v += (1 - c)*(Li3(r) + Li3(-r)) - 2*c*Li3(2*c/(1 + c))
    v += Li2((1 + c)/2)*((1 + c)*lp - (1 - c)*lm) + q*(Li2(c) - Li2(-c))*(lp + lm)
    v += -(1 - c)*lm**2*lp + (1 - c)*lm*lp**2 + (4*c**2 - c + 1)/(6*c)*lp**3 + q*l0*lm**2 - (1 - c)*l0*lm*lp \
        - (1 + c**2)/(2*c)*l0*lp**2
    v += lam*((1 - c)*lm**2 - (1 - c)*lm*lp - c*lp**2) + lam**2*(1 - c)/2*(lp - lm) \
        - mp.pi**2*((1 - c)/(6*c)*lm + (1 + c)/(6*c)*lp)
    v += (mp.mpf(7)/4 + 2*c)*z3 - c/6*mp.pi**2*lam + c/3*lam**3
    return v


# ---- own HPL implementation: power series in x, letters 0, 1, -1 (last letter nonzero)
def hpl_coeffs(word, N):
    """coefficients a[1..N] of H_word(x) = sum a_k x^k (word tuple, last letter nonzero)."""
    a = [mp.mpf(0)]*(N + 1)
    last = word[-1]
    for k in range(1, N + 1):
        a[k] = mp.mpf(1)/k if last == 1 else mp.mpf((-1)**(k + 1))/k
    for letter in reversed(word[:-1]):
        b = [mp.mpf(0)]*(N + 1)
        if letter == 0:
            for k in range(1, N + 1):
                b[k] = a[k]/k
        else:
            # f(t) H(t) = sum_m (sum_{k<=m} letter^(m-k) a_k) t^m, then integrate: t^m -> x^(m+1)/(m+1)
            acc = mp.mpf(0)
            for m in range(1, N):
                acc = acc*letter + a[m] if letter == -1 else acc + a[m]
                b[m + 1] = acc/(m + 1)
        a = b
    return a


def hpl(word, x, N):
    a = hpl_coeffs(word, N)
    s = mp.mpf(0); p = mp.mpf(1)
    for k in range(1, N + 1):
        p *= x
        s += a[k]*p
    return s


def S_hpl(x, N):
    H = lambda w: hpl(w, x, N)
    return (1 - x)*(H((0, 0, 1)) + H((0, 1, -1)) - H((1, -1, 1))) + (1 - x**2)/(2*x)*H((1, 0, 1)) \
        - (1 - x)**2/(2*x)*H((1, 0, -1))


def M2_hpl(c, N):
    q = (1 - c**2)/(2*c)
    return z2/2 - Li2(c**2)/2 - z2*q*mp.atanh(c) + S_hpl(c, N) + S_hpl(-c, N)


def W2_tint(a):
    # W2 = - int_1^{t1} log g(t) dt/(t+1), g = (t^2-1)/(2(t - cos a)), t1 = 1 + 2 sin(a/2)
    ca = mp.cos(a); t1 = 1 + 2*mp.sin(a/2)
    return -mp.quad(lambda t: mp.log((t**2 - 1)/(2*(t - ca)))/(t + 1), [1, (1 + t1)/2, t1])


def W3_x(a):
    s, c = mp.sin(a/2), mp.cos(a/2)
    tau = mp.tan(a/4)
    k = (c**2 + 1)/(4*c)
    table = [(0, (0, -1, mp.mpf(1)/2)), (1/s, (-1, 1, -mp.mpf(1)/2)), (s, (mp.mpf(1)/2, 0, mp.mpf(1)/2)),
             (tau, (-k, c, -c)), (1/tau, (k, -c, c))]

    def f(x):
        lg = mp.log(x/s) + mp.log(1 + s*x) - mp.log(1 + x/s)
        lm = mp.log(1 + x/tau) + mp.log(1 + tau*x) - mp.log(1 + x/s)
        tot = 0
        for l_, (A, B, C) in table:
            tot += (A*lg**2 + B*lg*lm + C*lm**2)/(x + l_)
        return tot
    return mp.quad(f, [0, mp.mpf(1)/8, mp.mpf(1)/2, 1])


def M2_2pi3():
    L2, L3 = mp.log(2), mp.log(3)
    return (Li2(mp.mpf(1)/3) + L2**2 - L2*L3 + L3**2/2 + mp.mpf(15)/2*Li3(mp.mpf(1)/3) + mp.mpf(17)/2*Li3(-mp.mpf(1)/2)
            + mp.mpf(3)/16*z3 - mp.mpf(5)/12*mp.pi**2*L2 + mp.pi**2*L3/4 - L2**3/12 + mp.mpf(5)/4*L2**2*L3
            - mp.mpf(7)/4*L2*L3**2 + mp.mpf(5)/8*L3**3 - mp.mpf(7)/2*Li2(mp.mpf(1)/3)*L2 + mp.mpf(15)/4*Li2(mp.mpf(1)/3)*L3)


if __name__ == '__main__':
    worst = {}
    for tag, a in ANG.items():
        V = rd(f'runB/vc_{tag}.txt')
        c = mp.cos(a/2)
        d = {}
        d['B0'] = abs(M0f(c) - V[0])
        d['B1'] = abs(M1f(c) - V[1])
        d['Phi'] = abs(Phi(c) - Phi_series(c))
        d['W2closed_vs_tint'] = abs(W2_tint(a) - (z2/2 - Li2(c**2)/2))
        d['W2+W3x'] = abs(z2/2 - Li2(c**2)/2 + W3_x(a) - V[2])
        d['B2prime'] = abs(M2_classical(c) - V[2])
        Nh = int(60*mp.log(10)/(-mp.log(c))) + 50
        d['B2_hpl'] = abs(M2_hpl(c, Nh) - V[2])
        print(f'{tag:>7} c={mp.nstr(c, 8)}: ' + '  '.join(f'{k}={mp.nstr(v, 3)}' for k, v in d.items()), flush=True)
        for k, v in d.items():
            worst[k] = max(worst.get(k, 0), v)
    print('WORST over 13 angles:', '  '.join(f'{k}={mp.nstr(v, 3)}' for k, v in worst.items()))
    V = rd('runB/vc_2pi_3.txt')
    print('2pi/3 special formula - M2 (Cauchy) =', mp.nstr(M2_2pi3() - V[2], 3))
    mp.mp.dps = 330
    z2, z3, lam = mp.zeta(2), mp.zeta(3), mp.log(2)   # recompute constants at 330 digits
    sol = rd('values_2pi_3.txt')
    print('2pi/3 special formula - solver 318-digit M2 =', mp.nstr(M2_2pi3() - sol[2], 3))
    c = mp.sqrt(3)/2
    print('pi/3: B2prime at c=sqrt3/2 (330 dps) - solver 318-digit M2 =', mp.nstr(M2_classical(c) - rd('values_pi_3.txt')[2], 3))
    for tag in ['pi_2', '3pi_4', 'pi_6', '5pi_6', 'pi_4', 'gen1p3', 'gen2p0']:
        c = mp.cos(ANG_hi[tag]()/2)
        print(f'{tag}: B2prime (330 dps) - solver 318-digit M2 =', mp.nstr(M2_classical(c) - rd(f'values_{tag}.txt')[2], 3),
              '; B1 - solver M1 =', mp.nstr(M1f(c) - rd(f'values_{tag}.txt')[1], 3))
