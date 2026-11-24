"""m2_split_check.py -- independent numerical checks of the exact split M2 = W2 + W3 (m2_structure.py).

x-variables (x = (t-1)/(2s) in [0,1], s = sin(a/2), c = cos(a/2), tau = tan(a/4)):
  lg = log g     = log(x/s) + log(1+s x) - log(1+x/s)
  lm = log(1+g)  = log(1+x/tau) + log(1+tau x) - log(1+x/s)
  W2 = - int_0^1 lg dx/(x+1/s)
  W3 = int_0^1 sum_sigma [ A_s lg^2 + B_s lg lm + C_s lm^2 ] dx/(x+lam_s)
  with (lam; A,B,C) = (0; 0,-1,1/2), (1/s; -1,1,-1/2), (s; 1/2,0,1/2),
                      (tau; -(c^2+1)/(4c), c, -c), (1/tau; (c^2+1)/(4c), -c, c)   [residues, m2_structure.out]
Checks: (i) W2 + W3 vs stored corner_series M2 (values_*.txt), (ii) W2 vs closed form
  W2 = 1/2 Li2(s^2) + (log terms), evaluated from the explicit antiderivative,
(iii) M2(pi/2) vs the 1200-digit squircle value, (iv) 3pi/4 and pi/6 vs D10b's 330-digit values.
"""
import re
import mpmath as mp

DPS = 130
mp.mp.dps = DPS + 20
ANG = {'pi_2': mp.pi/2, '2pi_3': 2*mp.pi/3, 'pi_3': mp.pi/3, '3pi_4': 3*mp.pi/4, 'pi_6': mp.pi/6,
       '5pi_6': 5*mp.pi/6, 'pi_4': mp.pi/4, 'gen1p3': mp.mpf(13)/10, 'gen2p0': mp.mpf(2)}


def stored(tag, k=2):
    for line in open(f'values_{tag}.txt'):
        if line.startswith(f'M{k} ='):
            return mp.mpf(line.split('=')[1].strip())


def parts(a):
    c, s = mp.cos(a/2), mp.sin(a/2)
    tau = mp.tan(a/4)
    q = (c*c + 1)/(4*c)
    rows = [(0, 0, -1, mp.mpf(1)/2), (1/s, -1, 1, -mp.mpf(1)/2), (s, mp.mpf(1)/2, 0, mp.mpf(1)/2),
            (tau, -q, c, -c), (1/tau, q, -c, c)]

    def lg(x):
        return mp.log(x/s) + mp.log1p(s*x) - mp.log1p(x/s)

    def lm(x):
        return mp.log1p(x/tau) + mp.log1p(tau*x) - mp.log1p(x/s)

    def f3(x):
        G, Mm = lg(x), lm(x)
        tot = 0
        for lam, A, B, C in rows:
            tot += (A*G*G + B*G*Mm + C*Mm*Mm)/(x + lam)
        return tot
    pts = [0, mp.mpf(1)/1024, mp.mpf(1)/32, mp.mpf(1)/4, 1]
    W3 = mp.quad(f3, pts)
    W2 = -mp.quad(lambda x: lg(x)/(x + 1/s), pts)
    return W2, W3


def W2_closed(a):
    """explicit antiderivative (t-variables): W2 = -I1 - Im1 + Ica + log2*log(1+s)."""
    c, s = mp.cos(a/2), mp.sin(a/2)
    I1 = mp.log(2*s)*mp.log(1 + s) + mp.polylog(2, -s)
    Im1 = (mp.log(2*(1 + s))**2 - mp.log(2)**2)/2
    Ica = (mp.log(2*s*(1 + s))*mp.log((1 + s)/c**2) + mp.polylog(2, -s*(1 + s)/c**2)
           - mp.log(2*s*s)*mp.log(1/c**2) - mp.polylog(2, -s*s/c**2))
    return -I1 - Im1 + Ica + mp.log(2)*mp.log(1 + s)


def W2_reduced(a):
    """after Landen: W2 = 1/2 Li2(s^2) + logs (see m2_w2_symbolic.py for the log part)."""
    c, s = mp.cos(a/2), mp.sin(a/2)
    l1p, l1m, ls, l2 = mp.log(1 + s), mp.log(1 - s), mp.log(s), mp.log(2)
    # direct transcription of Landen-reduced antiderivative (see docstring of m2_w2_symbolic.py)
    I1 = mp.log(2*s)*l1p + mp.polylog(2, -s)
    Im1 = (mp.log(2*(1 + s))**2 - l2**2)/2
    Ica = (mp.log(2*s*(1 + s))*mp.log((1 + s)/c**2) + (-mp.polylog(2, s) - l1m**2/2)
           - mp.log(2*s*s)*mp.log(1/c**2) - (-mp.polylog(2, s*s) - mp.log(1 - s*s)**2/2))
    return -I1 - Im1 + Ica + l2*l1p


if __name__ == '__main__':
    print(f'quadrature at dps {mp.mp.dps}')
    for tag, a in ANG.items():
        W2, W3 = parts(a)
        M2 = stored(tag)
        d1 = W2 + W3 - M2
        d2 = W2 - W2_closed(a)
        d3 = W2 - W2_reduced(a)
        print(f'{tag:>7}: W2 = {mp.nstr(W2, 25)}  W3 = {mp.nstr(W3, 25)}  |W2+W3-M2| = {mp.nstr(abs(d1), 3)}'
              f'  |W2-closed| = {mp.nstr(abs(d2), 3)}  |W2-reduced| = {mp.nstr(abs(d3), 3)}', flush=True)
    # (iii) squircle
    sq = mp.mpf(open('./paper/squircle_note/outputs/'
                     'M2_methodA_1200.txt').read().split()[0].strip())
    mp.mp.dps = 330
    print('pi/2 : |M2(series, 318 dig) - squircle M2 (1200 dig)| =',
          mp.nstr(abs(stored('pi_2') - sq), 3))
    txt = open('research/08_theory/D10b_octagon/M_new_340.txt').read()
    for tag, key in [('3pi_4', '3pi/4'), ('pi_6', 'pi/6')]:
        mm = re.search(r'M2\(' + re.escape(key) + r'\)\s*=\s*([-0-9.]+)', txt)
        if mm:
            v = mp.mpf(mm.group(1))
            print(f'{key}: |M2(ours) - M2(D10b, {len(mm.group(1))-3} dig)| =', mp.nstr(abs(stored(tag) - v), 3))
