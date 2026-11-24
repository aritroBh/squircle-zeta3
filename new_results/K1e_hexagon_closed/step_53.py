
import mpmath as mp, json, sys
WORK='research/08_theory/K1e_hexagon_closed'
mp.mp.dps = 50
sys.path.insert(0, 'research/10_frontier/E1_edge_law/verify')
from vjet import integrand_jet
import sympy as sp
# compare my arc-length integrand (d=1 value) with polar integrand * dtheta/dt
e = json.load(open(f"{WORK}/geometry_fresh.json"))[0]
A0=mp.mpf(sp.Rational(e['A0'])); A1=mp.mpf(sp.Rational(e['A1'])); B=mp.mpf(sp.Rational(e['B']))
kap=mp.mpf(sp.Rational(e['kap'])); L2=mp.mpf(sp.Integer(e['L2'])); q0=mp.mpf(sp.Integer(e['q0'])); V02=mp.mpf(sp.Integer(e['V02']))
w=[mp.mpf(sp.Rational(e['w0'])), mp.mpf(sp.Rational(e['w1'])), mp.mpf(sp.Rational(e['wo']))]
t=mp.mpf('0.5')
av=[A0*t, A1*(1-t), B*(1+kap*t)]
uv=[av[j]/w[j] for j in range(3)]
LA=[mp.log(u) for u in uv]
F0=mp.fsum([av[j]*uv[j] for j in range(3)])
F1=mp.fsum([av[j]*uv[j]*LA[j] for j in range(3)])
V0=sp.Matrix([1,0]); D=sp.Matrix([2,3])
yD=q0+t*L2; yy=V02+2*q0*t+L2*t*t
Phi0=F1/(2*F0)
G0 = F0**(-1)*(L2 - 2*yD*Phi0 + yy*Phi0**2)
print("my G0 at t=0.5:", mp.nstr(G0,15), " sqrt:", mp.nstr(mp.sqrt(G0),15))
# polar at the same point
vs = [(-mp.mpf(0)/18, mp.mpf(1)/18), (-mp.mpf(3)/18, mp.mpf(3)/18), (-mp.mpf(1)/18, mp.mpf(0)/18)]
vs = [(-mp.mpf(int(vv[1]))/18, mp.mpf(int(vv[0]))/18) for vv in [sp.Matrix([1,0]), sp.Matrix([3,3]), sp.Matrix([0,1])]]
cs=[mp.mpf(3),mp.mpf(5),mp.mpf(3)]; lams=[0,0,0]
yx = 1+2*t; yyv = 3*t
th = mp.atan2(yyv, yx)
pol = integrand_jet(th, vs, cs, lams)[0]
# dtheta/dt: th=atan2(3t,1+2t); dth/dt = (3(1+2t)-2*3t)/((1+2t)^2+(3t)^2) = 3/((1+2t)^2+9t^2)
dth = 3/((1+2*t)**2+9*t*t)
print("polar integrand*dth/dt:", mp.nstr(pol*dth,15))
print("ratio:", mp.nstr(mp.sqrt(G0)/(pol*dth),8))
