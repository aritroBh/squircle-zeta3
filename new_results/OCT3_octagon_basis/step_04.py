import json, time
import sympy as sp
from mpmath import mp, mpf, sqrt, log, pi, zeta, polylog, nstr

mp.dps = 280
t0=time.time()
r2s = sp.sqrt(2); bs = r2s-1
zs_sym = [ -1, -r2s/2, -bs, (1-r2s)/2, -bs**2, (4-3*r2s)/8, bs**4, (2-r2s)/4,
           bs**2, 1-r2s/2, bs, sp.Rational(1,2), 2-r2s, r2s/2, 2*bs, (2+r2s)/4, 12*r2s-16 ]
znames = ["-1","-r2/2","-b","(1-r2)/2","-b^2","(4-3r2)/8","b^4","(2-r2)/4",
          "b^2","1-r2/2","b","1/2","2-r2","r2/2","2b","(2+r2)/4","12r2-16"]
def to_mpf(e): return mpf(sp.N(e, 300))
zs = [to_mpf(z) for z in zs_sym]
for n,z in zip(znames, zs):
    assert -1 <= z < 1, (n,z)
print("all 17 z in [-1,1): OK   [%.1fs]"%(time.time()-t0))

Li3 = []
for n,z in zip(znames,zs):
    v = polylog(3, z)
    Li3.append(v)
    print("Li3(%-10s) = %s" % (n, nstr(v, 30)), " [%.1fs]"%(time.time()-t0))

r2 = sqrt(2); b = r2-1
l2 = log(2); lu = log(1+r2); z3 = zeta(3); pi2 = pi**2
base = list(Li3) + [z3, pi2*l2, pi2*lu, l2**3, l2**2*lu, l2*lu**2, lu**3,
                    pi2, l2**2, l2*lu, lu**2, l2, lu, mpf(1)]
basenames = ["Li3(%s)"%n for n in znames] + ["z3","pi2*l2","pi2*lu","l2^3","l2^2*lu","l2*lu^2","lu^3",
                    "pi2","l2^2","l2*lu","lu^2","l2","lu","1"]
assert len(base)==31

# extended x {1, sqrt2}
consts = []; names=[]
for nm,c in zip(basenames, base):
    consts += [c, r2*c]; names += [nm, "r2*"+nm]
print("total constants:", len(consts))

# sanity: known identity Li3(1/2) = 7/8 z3 - pi2/12 l2 + 1/6 l2^3
res = Li3[11] - (sp.Rational(7,8)*mpf(0)+ (7*z3/8 - pi2*l2/12 + l2**3/6))
print("sanity Li3(1/2) identity residual:", nstr(res, 5))
# sanity: Landen for b: Li3(b)+Li3(2-r2)+Li3(-r2/2) = z3 + pi2/6*log(1-b) - 1/2 log(b) log(1-b)^2 + 1/6 log(1-b)^3
lb = log(b); l1b = log(2-r2)
res2 = Li3[10]+Li3[12]+Li3[1] - (z3 + pi2*l1b/6 - lb*l1b**2/2 + l1b**3/6)
print("sanity Landen(b) residual:", nstr(res2, 5))
# sanity: duplication Li3(b^2) = 4(Li3(b)+Li3(-b))
res3 = Li3[8] - 4*(Li3[10]+Li3[2])
print("sanity Li3(b^2)=4(Li3(b)+Li3(-b)) residual:", nstr(res3, 5))

json.dump({"names":names, "vals":[nstr(c, 280) for c in consts]},
          open("consts62.json","w"))
json.dump({"basenames":basenames, "base":[nstr(c, 280) for c in base]},
          open("consts31.json","w"))
print("saved consts62.json / consts31.json   [%.1fs]"%(time.time()-t0))
