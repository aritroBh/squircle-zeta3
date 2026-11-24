import mpmath as mp
mp.mp.dps = 40

# brute force triple sum S = sum_{1<=n1<n2<n3<=M} i^{n3}/(n1 n2 n3)
def brute(M):
    s = mp.mpc(0)
    for n3 in range(3, M+1):
        inner = mp.mpf(0)
        for n2 in range(2, n3):
            for n1 in range(1, n2):
                inner += 1/(mp.mpf(n1)*n2)
        s += (1j)**n3 * inner / n3
    return s

for M in (100, 400):
    print(M, mp.nstr(brute(M), 20))
# (H^2-H2)/2 formula
def formula(M):
    s = mp.mpc(0); H = mp.mpf(0); H2 = mp.mpf(0)
    for n in range(1, M+1):
        s += (1j)**n * (H*H - H2)/2 / n
        H += 1/mp.mpf(n); H2 += 1/mp.mpf(n)**2
    return s
print('formula 400:', mp.nstr(formula(400), 20))
print('quad said: 0.0999539936646647 - 0.0335771478470153i')
