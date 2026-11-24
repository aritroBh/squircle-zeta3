import mpmath as mp
def _mul(A, B, N):
    return [mp.fsum(A[i]*B[n - i] for i in range(n + 1)) for n in range(N + 1)]
def _sqrt(Q, N):
    b = [mp.sqrt(Q[0])]
    for n in range(1, N + 1):
        s = Q[n] - mp.fsum(b[k]*b[n - k] for k in range(1, n))
        b.append(s/(2*b[0]))
    return b
def integrand_series(g, c, N):
    L = mp.log(g)
    m = mp.log1p(g)
    fac = [mp.factorial(k) for k in range(N + 1)]
    E1 = [(-L)**k/fac[k] for k in range(N + 1)]
    E2 = [(-2*L)**k/fac[k] for k in range(N + 1)]
    Q = [1 + g*g - 2*c*g] + [g*g*E2[k] - 2*c*g*E1[k] for k in range(1, N + 1)]
    B = _sqrt(Q, N)
    A = [(-m)**k/fac[k]/(1 + g) for k in range(N + 1)]
    AB = _mul(A, B, N)
    AB[0] -= 1
    G = [L**k/fac[k]/g for k in range(N + 1)]
    return _mul(AB, G, N)

stored = mp.mpf("0.010717362932289465418409576291224895855718975738464335760155612417841103619782418171396503742199126038436380620836337865304662784078999406937965739269009749615547401879802207222417678105864671523471037376024873555429679637766699958870900647844995948532283743592952407677723163070778394289744442351066844742744206712464681")
for dps in (50, 100, 200, 300):
    mp.mp.dps = dps
    c = mp.cos(mp.pi/3)
    cache = {}
    def F(g):
        if g not in cache:
            cache[g] = integrand_series(g, c, 3)
        return cache[g][3]
    v = mp.quad(F, [0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1])
    print(f"dps {dps}: M3 = {mp.nstr(v, 20)}  |err| = {mp.nstr(abs(v - stored), 3)}")
