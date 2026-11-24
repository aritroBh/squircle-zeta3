"""
H1: compare the exact Taylor coefficients (Route A/B json) with the direct quadrature values of S(p).
Usage: python 04_compare_series_vs_direct.py coeffs.json n direct_S_dpsXX.json
Prints, for each d: S_direct, series truncated at each order, residuals; and the
'peeled' remainders  R_j(d) = (S/S(1) - sum_{k<=j} c_k d^k)/d^{j+1}  next to c_{j+1}.
"""
import sys, json
import mpmath as mp
import sympy as sp

mp.mp.dps = 60
cfile, n, dfile = sys.argv[1], int(sys.argv[2]), sys.argv[3]
C = json.load(open(cfile))
C = C[str(n)] if str(n) in C else C
dcoef_str = C['d']
K = len(dcoef_str) - 1
zs = {sp.Symbol(f'z{j}'): sp.zeta(j) for j in range(2, 40)}
c = [mp.mpf(sp.N(sp.sympify(s).subs(zs), 70)) for s in dcoef_str]
S1 = 2 ** n * mp.sqrt(n) / mp.factorial(n - 1)
D = json.load(open(dfile))
print(f'n={n}, K={K}, S_n(1)={mp.nstr(S1, 30)}')
print('coefficients c_k of S(1+d)/S(1):')
for k, ck in enumerate(c):
    print(f'  c_{k} = {mp.nstr(ck, 40)}')
print()
print(f'{"d":>8} {"S_direct":>44} {"S_series(K)":>44} {"diff":>10} {"|c_(K+1)| d^(K+1) scale":>12}')
for dstr, v in D.items():
    if dstr in ('0', '1'):
        continue
    d = mp.mpf(dstr)
    Sd = mp.mpf(v['N1'])
    ser = S1 * sum(ck * d ** k for k, ck in enumerate(c))
    print(f'{dstr:>8} {mp.nstr(Sd, 40):>44} {mp.nstr(ser, 40):>44} {mp.nstr(Sd - ser, 3):>10} {mp.nstr(abs(d) ** (K + 1), 2):>12}')
print()
print('peeled remainders R_j(d) = (S/S(1) - sum_{k<=j} c_k d^k)/d^{j+1}, which should tend to c_{j+1}:')
for j in range(0, 5):
    row = []
    for dstr in ('-0.001', '0.001'):
        d = mp.mpf(dstr)
        r = (mp.mpf(D[dstr]['N1']) / S1 - sum(c[k] * d ** k for k in range(j + 1))) / d ** (j + 1)
        row.append(mp.nstr(r, 15))
    # symmetric average removes the c_{j+2} d term
    d = mp.mpf('0.001')
    rp = (mp.mpf(D['0.001']['N1']) / S1 - sum(c[k] * d ** k for k in range(j + 1))) / d ** (j + 1)
    rm = (mp.mpf(D['-0.001']['N1']) / S1 - sum(c[k] * (-d) ** k for k in range(j + 1))) / (-d) ** (j + 1)
    print(f'  j={j}: R_j(-1e-3)={row[0]}, R_j(+1e-3)={row[1]}, avg={mp.nstr((rp + rm) / 2, 15)}, '
          f'c_{j + 1}={mp.nstr(c[j + 1], 15)}, avg-c_(j+1)={mp.nstr((rp + rm) / 2 - c[j + 1], 3)} (expect ~c_(j+3)*1e-6)')
