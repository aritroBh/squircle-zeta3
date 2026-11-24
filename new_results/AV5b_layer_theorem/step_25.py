
import mpmath as mp
mp.mp.dps = 60
def P_squircle(p):
    q = p/(p-1)
    return 8*mp.quad(lambda th: (mp.cos(th)**q + mp.sin(th)**q)**(1/q), [0, mp.pi/4])
for d in [mp.mpf('0.02'), mp.mpf('0.03')]:
    exact = P_squircle(1+d) - 4*mp.sqrt(2)
    print(f"d={mp.nstr(d,3)}: exact squircle DeltaP = {mp.nstr(exact, 35)}")
