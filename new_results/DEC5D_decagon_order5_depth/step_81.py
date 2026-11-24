
import mpmath as mp
mp.mp.dps = 250
tests = [((5,5,0), -mp.mpf(120)/6**6), ((0,0,2), mp.mpf(2)), ((4,0,5), -mp.mpf(120)/5**6), ((0,5,0), -mp.mpf(120)/1**6)]
for (n,a,b), ex in tests:
    f = lambda x: x**n * mp.log(x)**a * mp.log(1-x)**b
    v = mp.quad(f, [0,1], method='tanh-sinh')
    print((n,a,b), "err:", mp.nstr(abs(v-ex), 3))
# so quadrature is the issue. try higher maxdegree
f = lambda x: x**5 * mp.log(x)**5
v1 = mp.quad(f, [0,1], method='tanh-sinh')
mp.mp.dps = 300
v2 = mp.quad(f, [0,1], method='tanh-sinh')
ex = -mp.mpf(120)/6**6
mp.mp.dps = 250
print("deg default err(250dps):", mp.nstr(abs(v1-ex),3), " 300dps:", mp.nstr(abs(v2-ex),3))
