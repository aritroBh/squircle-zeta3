import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from mpmath import mp, mpf, pi, sqrt, exp, log
mp.dps = 60
import dcy2_audit as A

N = 400
tq = A.true_quant(N)
It = tq['It']   # direct ray quadrature (2/pi)*quad
print("It (direct)   =", mp.nstr(It, 12))

# window decomposition
y=tq['y']; Delta=tq['Delta']; rho=1/y
s0 = tq['s0']
phi_s0 = (N-1)*log(s0)-N*log(1+s0)+1j*pi*s0/2
ellt = tq['ellt']; F = tq['F']
W = mp.quad(lambda t: exp(ellt(t))*F(t), [-Delta, 0, Delta])
TR = mp.quad(lambda t: exp(ellt(t))*F(t), [Delta, rho, 10*rho])
TL = mp.quad(lambda t: exp(ellt(t))*F(t), [-rho+mpf('1e-9'), -Delta])
Jw = (2/pi)*exp(phi_s0)*A.omega*(W+TR+TL)
print("J (window)    =", mp.nstr(Jw, 12))
print("diff          =", mp.nstr(abs(It-Jw), 6))
print("|It| =", mp.nstr(abs(It),6), " e^{-theta}=", mp.nstr(exp(-sqrt(pi*N)),6))
# pieces
print("|W| =", mp.nstr(abs(W),6), " |TR| =", mp.nstr(abs(TR),6), " |TL| =", mp.nstr(abs(TL),6))
