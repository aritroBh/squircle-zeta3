code = r'''import mpmath as mp
mp.mp.dps = 50
w = mp.expjpi(mp.mpf(1)/4)
a = mp.pi/(2*mp.sqrt(2))
def H(s): return mp.sqrt(mp.pi)/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf(3)/2)
def K(s): return 1/(2*mp.cos(mp.pi*s/2))
def E(s): return mp.expj(mp.pi*s/2)
print("CX6R step_04b: residues + (13), 50 dps")
for n in range(4):
    s0 = mp.mpf(2*n+1)
    res = mp.limit(lambda s: (s-s0)*K(s), s0)  # simple-pole residue
    print(f"Res K at {s0}: {mp.nstr(res,20)}  vs (-1)^(n+1)/pi = {mp.nstr((-1)**(n+1)/mp.pi,20)}")
m1 = m2 = m3 = mp.mpf(0)
for i in range(4000):
    r = mp.mpf('0.005')*i
    s = mp.mpf('0.5')+r*w
    m1 = max(m1, abs(H(s)))
    m2 = max(m2, abs(K(s))*mp.exp(a*r))
    m3 = max(m3, abs(K(s)-E(s))*mp.exp(3*a*r))
print("max |H| on Gamma grid (bound 2):", mp.nstr(m1,10))
print("max |K| e^{ar} (bound 2):", mp.nstr(m2,10))
print("max |K-E| e^{3ar} (bound 2):", mp.nstr(m3,10))
mn = min(abs(1+mp.exp(mp.j*mp.pi*(mp.mpf('0.5')+mp.mpf('0.0005')*i*w))) for i in range(40000))
print("min |1+e^{i pi s}| on Gamma (proof uses >= 1/2):", mp.nstr(mn,15))
# also: |K(2J+iy)| = 1/(2 cosh(pi y/2)) check (vertical segment, sec.2)
d = max(abs(abs(K(2*J+mp.j*y)) - 1/(2*mp.cosh(mp.pi*y/2)))
        for J in (1,3,7) for y in (mp.mpf('0.3'), mp.mpf('2.5'), mp.mpf('11')))
print("vertical-segment |K| formula max err:", mp.nstr(d,6))
# int_0^inf sech(pi y/2) dy = 1 (used for the O(1/J) bound)
print("int sech(pi y/2) dy over [0,inf):", mp.nstr(mp.quad(lambda y: 1/mp.cosh(mp.pi*y/2), [0,20,mp.inf]),20))
'''
open('step_04b.py','w').write(code)
exec(compile(code,'step_04b.py','exec'))
