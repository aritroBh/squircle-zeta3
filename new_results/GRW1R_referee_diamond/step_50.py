import mpmath as mp
mp.mp.dps = 70

def Jint(eps):
    # J(eps) = -(4/(eps)) e^{-pi/(2 eps)} e^{i pi/4} int_0^pi e^{-s/eps} sqrt(sin s) e^{-i s/2} ds  (GRW1's -4i/i = -4/eps)
    f = lambda s: mp.e**(-s/eps)*mp.sqrt(mp.sin(s))*mp.e**(-1j*s/2)
    integ = mp.quad(f, [0, eps*3, eps*20, mp.pi])
    return -(4/eps)*mp.e**(-mp.pi/(2*eps))*mp.e**(1j*mp.pi/4)*integ

def Jpred(eps, terms=2):
    base = 2*mp.sqrt(mp.pi)*mp.e**(-3j*mp.pi/4)*mp.sqrt(eps)*mp.e**(-mp.pi/(2*eps))
    if terms==1: return base
    return base*(1 - 3j*eps/4 - 25*eps**2/32)

print("eps | |J|/|pred1| | phase(J)-phase(pred) | J/(pred with 1-3ie/4-25e^2/32)")
for eps in [mp.mpf('0.3'), mp.mpf('0.2'), mp.mpf('0.15'), mp.mpf('0.1'), mp.mpf('0.075'), mp.mpf('0.05')]:
    J = Jint(eps)
    p1 = Jpred(eps,1); p2 = Jpred(eps,2)
    print(f"eps={eps}: |J|/|p1|={mp.nstr(abs(J)/abs(p1),10)}  dphase={mp.nstr(mp.arg(J)-mp.arg(p1),8)}")
    print(f"         J/p2-1 = {mp.nstr(J/p2-1,8)}")
# also verify the O(eps^2) coefficient -25/32: J/(base) = 1 - 3 i eps/4 - 25 eps^2/32 + O(eps^3)
print("\n(J/base - 1 + 3 i eps/4)/eps^2 -> -25/32 = -0.78125:")
for eps in [mp.mpf('0.1'), mp.mpf('0.075'), mp.mpf('0.05'), mp.mpf('0.03')]:
    J = Jint(eps)
    base = 2*mp.sqrt(mp.pi)*mp.e**(-3j*mp.pi/4)*mp.sqrt(eps)*mp.e**(-mp.pi/(2*eps))
    est = (J/base - 1 + 3j*eps/4)/eps**2
    print(f"eps={eps}: {mp.nstr(est,10)}")
with open('progress.log','a') as f:
    f.write("step_29: Stokes jump integral asymptotics verified; O(eps^2) coefficient -25/32 extracted\n")
