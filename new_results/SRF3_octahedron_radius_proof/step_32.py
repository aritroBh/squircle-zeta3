
import os, math
import mpmath as mp
mp.mp.dps = 50
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
print("engine log tail:", open(wd+"/engine_K120.log").read()[-200:])

# Fitted values (K=80, step above): C=1.376611, phi=-2.24764938
C_fit = mp.mpf("1.376611")
phi_fit = mp.mpf("-2.24764938")
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = mp.arg(z0)
# Theory: D = (sqrt(3)*pi*i/2) * e^{z0} * z0^{1/2}  [Borel singular part D*(1-z/z0)^{1/2} for S/S(1)]
D_pred = (mp.sqrt(3)*mp.pi*1j/2) * mp.exp(z0) * z0**mp.mpf('0.5')
C_pred = abs(D_pred)/mp.sqrt(mp.pi)
print("\n|D_pred| =", mp.nstr(abs(D_pred),15), "  C_pred = |D|/sqrt(pi) =", mp.nstr(C_pred,15))
print("C_fit  =", mp.nstr(C_fit,15), " ratio =", mp.nstr(C_fit/C_pred,10))
# phase: c_k ~ -(1/sqrt(pi)) k^{-3/2} Re[D z0^{-k}] ; compare to C cos(th k + phi)
# Re[D e^{-ik th}] = |D| cos(k th - arg D);  so c_k ~ -(1/sqrt(pi)) k^{-3/2} |D| cos(k th - arg D)
# fitted: +C k^{-3/2} cos(th k + phi):  phi = -arg D + pi  (mod 2pi)
phi_pred = (-mp.arg(D_pred) + mp.pi) % (2*mp.pi)
if phi_pred > mp.pi: phi_pred -= 2*mp.pi
print("phi_fit =", mp.nstr(phi_fit,10), " phi_pred =", mp.nstr(phi_pred,10), " diff =", mp.nstr(phi_fit-phi_pred,6))
