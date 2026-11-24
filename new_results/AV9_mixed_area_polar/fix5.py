
import mpmath as mp
mp.mp.dps = 70
def Vnum(etaval):
    q = 1+mp.mpf(etaval)
    if q < 3:   # v-form is fine
        f = lambda v: (2/q)*(v*(1-v))**(1/q-1)*(v**(q-1)+(1-v)**(q-1))**(1/q)
        return mp.quad(f, [0, mp.mpf('0.5'), 1])
    else:       # theta-form well conditioned
        p = q/(q-1)
        f = lambda th: 2*(q-1)*(mp.cos(th)**p+mp.sin(th)**p)**(1/p)*(mp.cos(th)*mp.sin(th))**(q-2)*(mp.cos(th)**q+mp.sin(th)**q)**(1/q-2)
        return mp.quad(f, [0, mp.pi/4, mp.pi/2])
print("Duality V(eta) vs V(1/eta), well-conditioned:")
for ev in ['0.25','0.1','0.05','1','3','0.5']:
    a = Vnum(ev); b = Vnum(mp.mpf(1)/mp.mpf(ev))
    print(f"eta={ev}: V={mp.nstr(a,25)}  V(1/eta)={mp.nstr(b,25)}  diff={mp.nstr(abs(a-b),3)}")
