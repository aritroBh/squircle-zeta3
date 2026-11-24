# Tasks (2 cont.) and (3): assemble E_{5,m} from MY e's (my_e5m.json) and MY N's (n_values110.json),
# compare with e5m_values.json / n6_values.json; E_{5,5}/24 vs c(6); LEAD law.
import json, sys
from fractions import Fraction
import mpmath as mp

BASE = '<project>/research/08_theory/WT6S_weight6_symbol_all_k/'
DPS = 110
mp.mp.dps = DPS + 25
s2 = mp.sqrt(2)

def fr(x):
    x = Fraction(x)
    return mp.mpf(x.numerator) / x.denominator

myE = json.load(open('my_e5m.json'))
Nj = json.load(open('n_values110.json'))
NA = {k: mp.mpf(v) for k, v in Nj['A'].items()}
NB = {k: mp.mpf(v) for k, v in Nj['B'].items()}

# compare to n6_values.json (the job's own N)
n6 = json.load(open(BASE + 'n6_values.json'))
print('n6_values.json keys:', list(n6.keys())[:6], '...')

# my N (A and B) vs the job's n6_values.json (its keys '5,0|0' <-> mine '50|0')
worstN = 0
for k, v in NA.items():
    kk = k[0] + ',' + k[1:]
    d = abs(v - mp.mpf(n6['values'][kk]))
    worstN = max(worstN, d)
print('worst |N_mine(A) - n6_values.json| over 12 integrals:', mp.nstr(worstN, 4))
print('worst |N_mine(A) - N_mine(B)| (re-read from json):', mp.nstr(max(abs(NA[k]-NB[k]) for k in NA), 4))

cells = list(myE.keys())            # '50|0', ...
def Evals(Ntab):
    out = {}
    for m in range(1, 6):
        tot = mp.mpf(0)
        for cell in cells:
            key, pole = cell.split('|')
            e = fr(myE[cell]['rat'][m - 1]) + s2 * fr(myE[cell]['sqrt2'][m - 1])
            tot += e * Ntab[f'{key}|{pole}']
        out[m] = tot
    return out

EA = Evals(NA)
EB = Evals(NB)
ref = json.load(open(BASE + 'e5m_values.json'))['values']
print('--- E_{5,m} from my quadrature (A) vs e5m_values.json ---')
for m in range(1, 6):
    r = mp.mpf(ref[f'E5,{m}'])
    print(f'E5,{m} mine(A) =', mp.nstr(EA[m], 60))
    print(f'      mine(B) - mine(A) = {mp.nstr(EB[m]-EA[m], 3)}   |mine(A) - e5m_values.json| = {mp.nstr(abs(EA[m]-r), 3)}')

# Compare my N values with the job's n6_values.json, whatever format it has
try:
    sample = n6[list(n6.keys())[0]]
    print('n6_values.json first entry type:', type(sample))
except Exception as ex:
    print('n6 err', ex)

# ---- c(6) and the LEAD law ----
c_task = mp.mpf('0.33568755822402035700767227984351182721')
c_ref = mp.mpf(json.load(open(BASE + 'c6_lead_reference.json'))['value'])
c_mine = EA[5] / 24
print('--- E_{5,5}/24 vs c(6) ---')
print('E5,5/24 mine   =', mp.nstr(c_mine, 70))
print('task c(6) (38d) =', '0.33568755822402035700767227984351182721')
print('|mine - task c(6)| =', mp.nstr(abs(c_mine - c_task), 3), ' (task value has 38 sig digits, rounding bound ~5e-39)')
print('|mine - c6_lead_reference.json (405d)| =', mp.nstr(abs(c_mine - c_ref), 3))

# LEAD law: c(6) = sigma4 * A5 / (2*4!) * [y^6] E(y),  E(y)=Gamma(1+y/2)Gamma(1/2-y/2)/sqrt(pi)
# Series of log E via zeta values / polygamma at 1/2, then exp.
N = 8
l = [mp.mpf(0)] * (N + 1)
# log Gamma(1+z) = -gamma z + sum_{n>=2} (-1)^n zeta(n) z^n / n,  z = y/2
for n in range(1, N + 1):
    coef = (-mp.euler if n == 1 else (-1) ** n * mp.zeta(n) / n) * mp.mpf(2) ** (-n)
    l[n] += coef
# log Gamma(1/2 + w) - log Gamma(1/2) = sum_{n>=1} psi^{(n-1)}(1/2) w^n / n!,  w = -y/2
for n in range(1, N + 1):
    l[n] += mp.polygamma(n - 1, mp.mpf(1) / 2) / mp.factorial(n) * (-mp.mpf(1) / 2) ** n
e = [mp.mpf(0)] * (N + 1)
e[0] = mp.mpf(1)
for n in range(1, N + 1):
    e[n] = sum(kk * l[kk] * e[n - kk] for kk in range(1, n + 1)) / n
print('E(y) coefficients [y^0..y^8]:', [mp.nstr(x, 12) for x in e])
# cross-check with numerical Taylor at lower precision
with mp.workdps(60):
    Ef = lambda y: mp.gamma(1 + y / 2) * mp.gamma(mp.mpf(1) / 2 - y / 2) / mp.sqrt(mp.pi)
    tay = mp.taylor(Ef, 0, 6)
    print('mp.taylor [y^6] =', mp.nstr(tay[6], 30), '  series [y^6] =', mp.nstr(e[6], 30),
          '  diff =', mp.nstr(abs(tay[6] - e[6]), 3))
sigma4, A5 = 1, 16
c_lead = sigma4 * A5 / mp.mpf(2 * 24) * e[6]
print('LEAD c(6) = sigma4*A5/(2*4!)*[y^6]E =', mp.nstr(c_lead, 70))
print('|LEAD law - E5,5/24 (mine)| =', mp.nstr(abs(c_lead - c_mine), 3))
print('|LEAD law - c6_lead_reference.json| =', mp.nstr(abs(c_lead - c_ref), 3))
print('|LEAD law - task c(6)| =', mp.nstr(abs(c_lead - c_task), 3))
print('E5,5 itself = 24*c =', mp.nstr(EA[5], 60), '  vs 24*LEAD =', mp.nstr(24 * c_lead, 60))
# exact digit-agreement count
d = abs(c_mine - c_lead)
print('agreement digits (mine vs LEAD law): ~', int(-mp.log10(d)) if d != 0 else 'inf')
json.dump({'E_A': {str(m): mp.nstr(EA[m], 120) for m in EA}, 'E_B': {str(m): mp.nstr(EB[m], 120) for m in EB},
           'c_lead': mp.nstr(c_lead, 120)}, open('E_values_mine.json', 'w'), indent=1)
