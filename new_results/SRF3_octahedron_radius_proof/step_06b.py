
"""SRF3 step_06b: fiber min via continuation from w=1/2; K-factor bounds."""
import mpmath as mp
mp.mp.dps = 50
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0)
print("== (2) corner fiber min |z(w)|, e^{-2zw}+e^{-2z(1-w)}=-1, by continuation from w=1/2")
# at w=1/2: z = log(2) + i pi  (and conjugates). Continue in w.
for sgn in [1]:
    w = mp.mpf('0.5')
    z = mp.log(2)+sgn*1j*mp.pi
    for target in [mp.mpf('0.4'),mp.mpf('0.3'),mp.mpf('0.2'),mp.mpf('0.1'),mp.mpf('0.05'),mp.mpf('0.02'),mp.mpf('0.005')]:
        dw = (target-w)/20
        for _ in range(20):
            w += dw
            for _ in range(60):
                F = mp.e**(-2*z*w)+mp.e**(-2*z*(1-w))+1
                dF = -2*w*mp.e**(-2*z*w)-2*(1-w)*mp.e**(-2*z*(1-w))
                z = z-F/dF
        # now at w=target
        print(f"  w={mp.nstr(w,6)}: |z|={mp.nstr(abs(z),10)}  |z|/|z0|={mp.nstr(abs(z)/R,8)}  z={mp.nstr(z,12)}")
    print("  (analytic limits: w=1/2 -> 2|z0| =", mp.nstr(2*R,10), "; w->0 -> |z0| =", mp.nstr(R,10), ")")

print("\n== (3) K-factor bound |k_j(V)| vs e^{-jV}/j! : (1-e^{-V})^{-2 t} = sum k_j t^j")
def kcoef(V,j):
    L = -2*mp.log(1-mp.e**(-V))
    return L**j/mp.factorial(j)
for V in [mp.mpf('0.5'), mp.mpf('1'), mp.mpf('2'), mp.mpf('4')]:
    worst = max(abs(kcoef(V,j))*mp.factorial(j)/mp.e**(-j*V) for j in range(1,9))
    print(f"  V={V}: max_j |k_j| j! / e^{{-jV}} = {mp.nstr(worst,6)}  (<=1 required, need |log(1-e^-V)|<=e^-V: {mp.log(1-mp.e**(-V))<mp.e**(-V)})")
# the needed inequality: -log(1-e^{-V}) <= e^{-V}? at V=0.5: -log(0.3935)=0.9328 vs e^{-0.5}=0.6065: FAILS!
# so the clean R2 bound |k_j|<=e^{-jV}/j! is false here: we need (2|log(1-e^{-V})|)^j/j! <= (C e^{-V})^j/j!
# i.e. 2|log(1-e^{-V})| <= C e^{-V}: C(V)=2|log(1-e^{-V})|e^{V}:
for V in [mp.mpf('0.2'),mp.mpf('0.5'), mp.mpf('1'), mp.mpf('2'), mp.mpf('4'), mp.mpf('8')]:
    print(f"  V={V}: C(V)=2|log(1-e^{{-V}})|e^V = {mp.nstr(2*abs(mp.log(1-mp.e**(-V)))*mp.e**V,8)}")
