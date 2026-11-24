import mpmath as mp, pickle, sys
sys.path.insert(0,'.')
import ftfield as F
mp.mp.dps = 80
m = 7
g = pickle.load(open('geom14.pkl','rb'))
LL = mp.sqrt(F.feval(g['L'],80))
def fsum(th, p):
    return sum(mp.exp(p*mp.log(abs(mp.cos(th - i*mp.pi/m)))) for i in range(m))
def r_and_rp(th, p):
    f = fsum(th,p)
    fp = sum(-mp.sin(th - i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    return f**(-1/p), -(1/p)*f**(-1/p-1)*fp
def integrand(th, p):
    r, rp = r_and_rp(th,p); return mp.sqrt(r*r+rp*rp)
# periodicity check of integrand
for ths in ['0.03','0.11']:
    th=mp.mpf(ths)
    print("period check th=%s:"%ths, mp.nstr(integrand(th,1)-integrand(th+mp.pi/7,1),3))
# gold chamber (offset, 1e-9 sliver removed) vs centered chamber, at p=1
pc_gold = 2*m*mp.quad(lambda th: integrand(th,1), [mp.pi/2+mp.mpf('1e-9'), mp.pi/2+mp.pi/m-mp.mpf('1e-9')])
pc_cent = 2*m*mp.quad(lambda th: integrand(th,1), [-mp.pi/14, mp.pi/14])
print("P(1) gold chamber   =", mp.nstr(pc_gold,30), " diff vs 14L:", mp.nstr(pc_gold-14*LL,3))
print("P(1) centered       =", mp.nstr(pc_cent,30), " diff vs 14L:", mp.nstr(pc_cent-14*LL,3))
# integrand near the vertex theta=pi/2 (gold endpoint) vs near 0 (centered endpoint)
print("integrand at pi/2 (p=1):", mp.nstr(integrand(mp.pi/2,1),20))
print("integrand at 0    (p=1):", mp.nstr(integrand(mp.mpf('0'),1),20))
print("integrand at pi/2+pi/7 (p=1):", mp.nstr(integrand(mp.pi/2+mp.pi/7,1),20))
# behavior of integrand in p near the vertex (theta=pi/2): is it stiff?
for ps in ['1.0','1.001','0.999']:
    print("integrand(pi/2+1e-9, p=%s) = %s"%(ps, mp.nstr(integrand(mp.pi/2+mp.mpf('1e-9'),mp.mpf(ps)),20)))
