import mpmath as mp
mp.mp.dps = 105
def I_quad(m, j):
    return mp.quad(lambda s: s**m * mp.log1p(mp.e**(-s))**j, [0, 0.5, 2, 8, 32, mp.inf])
def I_stirling(m, j):
    # log^j(1+x) = j! sum_l s(l,j) x^l / l! (SIGNED Stirling I); I = j! m! sum s(l,j)/(l! l^{m+1})
    tot = mp.mpf(0)
    for l in range(j, 400):
        s = mp.stirling1(l, j)   # signed
        tot += mp.mpf(int(s)) / (mp.factorial(l) * mp.mpf(l)**(m+1))
    return mp.factorial(j)*mp.factorial(m)*tot
def eta(s): return (1-mp.power(2,1-s))*mp.zeta(s)
print(f"{'(m,j)':8} {'quadrature':30} {'Stirling-series':30} {'diff':10}")
for (m,j) in [(0,1),(1,1),(2,1),(3,1),(4,1),(0,2),(1,2),(2,2),(3,2),(0,3),(1,3),(2,3),(0,4)]:
    a, b = I_quad(m,j), I_stirling(m,j)
    print(f"({m},{j})    {mp.nstr(a,25)}  {mp.nstr(b,25)}  {mp.nstr(abs(a-b),2)}")
print("\nclosed forms: I_{m,1} = m! eta(m+2):")
for m in range(0, 6):
    print(f"  m={m}: {mp.nstr(I_quad(m,1) - mp.factorial(m)*eta(m+2), 3)}")
print("I_{0,2} = zeta(3)/4:", mp.nstr(I_quad(0,2) - mp.zeta(3)/4, 3))
# A_{2m} classical identity (needed for odd tops): 2 A_{2m} = (2m+1)eta(2m+1) - zeta(2m+1) - 2 sum_{j=1}^{m-1} eta(2j) zeta(2m+1-2j)
print("\nA_{2m} identity, A_s = sum (-1)^{l-1} H_l/l^s:")
def A_alt(s):
    return mp.nsum(lambda l: (-1)**(l-1)*mp.harmonic(l)/mp.mpf(l)**s, [1, mp.inf], tol=mp.mpf('1e-95'))
for m in range(1, 5):
    s = 2*m
    cf = ((s+1)*eta(s+1) - mp.zeta(s+1))/2 - mp.fsum(eta(2*j)*mp.zeta(s+1-2*j) for j in range(1, m))
    A = A_alt(s)
    print(f"  m={m}: |A_{s} - closed| = {mp.nstr(abs(A-cf),3)}")
import json
Itab = {f"{m},{j}": str(I_quad(m,j)) for m in range(0, 9) for j in range(1, 6)}
json.dump(Itab, open('I_mj_quad.json','w'))
