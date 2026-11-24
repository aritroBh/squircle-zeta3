
import mpmath as mp
mp.mp.dps = 200
import json
def I_u(eta, pts):
    f2 = lambda u: mp.exp(-u)*(2*mp.exp(-eta*u/2)*mp.sqrt(mp.cosh(eta*u)))*(1+mp.exp(-u))**(-2-eta)
    return mp.quad(f2, pts)
refs = {}
for eta in ['0.1','0.2','0.3','-0.1','-0.2','-0.3']:
    e = mp.mpf(eta)
    v1 = I_u(e, [0, 1, 5, 20, 100])
    v2 = I_u(e, [0, 2, 9, 30, 150])
    print(f"eta={eta}: agree {mp.nstr(abs(v1-v2),3)}")
    refs[eta] = mp.nstr(v1, 150)
    print("  ", mp.nstr(v1, 90))
json.dump(refs, open("iref200.json","w"))
