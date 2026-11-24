import mpmath as mp
mp.mp.dps = 80
def eta(s): return (1-mp.power(2,1-s))*mp.zeta(s)
S1 = mp.nsum(lambda l: (-1)**l*mp.harmonic(l-1)/mp.mpf(l)**4, [2, mp.inf], tol=mp.mpf('1e-70'))
S2 = mp.nsum(lambda l: (-1)**l*(mp.harmonic(l)-1/l)/mp.mpf(l)**4, [1, mp.inf], tol=mp.mpf('1e-70'))
A4 = mp.nsum(lambda l: (-1)**(l-1)*mp.harmonic(l)/mp.mpf(l)**4, [1, mp.inf], tol=mp.mpf('1e-70'))
print("4*sum (-1)^l H_{l-1}/l^4 =", mp.nstr(4*S1, 30))
print("4*(eta(5)-A_4)         =", mp.nstr(4*(eta(5)-A4), 30))
print("I_{2,2} quadrature      = 0.195745588199876253442971")
print("eta(5)-A4 - S1 =", mp.nstr(eta(5)-A4-S1, 10))
# the correct expansion of I_{m,2}: 2 m! sum s(l,2)/(l! l^{m+1}); m=2 -> l^{3}, and (l-1)!/l!=1/l -> 1/l^4. So both should agree...
# unless I_{m,j} formula has l^{m+1} with m=2 -> l^3, i.e. weight m+j+1 = 5: sum s(l,2)/(l! l^{m+1}): weight = (m+1)+(j)... 
# check K3's odd top for n=2 (D~_{5,1}) with the QUADRATURE values:
e13, e12 = mp.mpf('0.5'), mp.mpf('1')
D51 = 2*(e13*I31 + e12*I22/2) if False else None
I31 = mp.quad(lambda s: s**3*mp.log1p(mp.e**(-s)), [0,1,4,16,mp.inf])
I22 = mp.quad(lambda s: s**2*mp.log1p(mp.e**(-s))**2, [0,1,4,16,mp.inf])
D51 = 2*(mp.mpf('0.5')*I31 + mp.mpf('1')*I22/mp.factorial(2))
print("D~_{5,1} from quadrature =", mp.nstr(D51, 30))
print("K3 closed form 2!*[z5 + 2 eta(2) z3] =", mp.nstr(2*(mp.zeta(5) + 2*eta(2)*mp.zeta(3)), 30))
