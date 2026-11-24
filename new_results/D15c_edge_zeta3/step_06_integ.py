
import sympy as sp
import sym_G, time
t0=time.time()
G21, G12, syms = sym_G.build()
s0,s1,c0,c1,S0,Sc = syms['s0'],syms['s1'],syms['c0'],syms['c1'],syms['S0'],syms['Sc']
sf,cf,Sh,Ch = sp.symbols('sf cf Sh Ch')
subsd = {
  sp.sqrt(1+Sc**2/S0**2): 1/cf,
  s0: sf*Ch+cf*Sh, c0: cf*Ch-sf*Sh,
  s1: Sh*cf-Ch*sf, c1: -(Ch*cf+Sh*sf),
  S0: cf/Sh, Sc: -sf/Sh,
}
def prep(G):
    Ge = sp.expand(G)
    # substitute the sqrt first (it appears as a sub-expression power)
    Ge = Ge.subs(subsd)
    return sp.expand(Ge)
g21 = prep(G21); g12 = prep(G12)
gsum = sp.expand(g21+g12)
# Use sf^2+cf^2=1 to reduce: replace sf^2 -> 1-cf^2 repeatedly? Keep both for now.
num, den = sp.fraction(sp.together(gsum))
num = sp.expand(num); den = sp.expand(den)
print("time prep", time.time()-t0)
print("DEN =", sp.factor(den))
# Express num as poly in sf,cf,Sh,Ch
P = sp.Poly(num, sf,cf,Sh,Ch)
print("num monomial count:", len(P.monoms()))
deg_sf = max(mm[0] for mm in P.monoms())
deg_cf = max(mm[1] for mm in P.monoms())
print("max deg sf", deg_sf, "max deg cf", deg_cf)
# Since integrand even in phi (sf->-sf symmetric), keep only even sf powers
import sympy
x = sp.Symbol('x')
# Build integrand as function of phi and integrate each monomial
phi = sp.Symbol('phi')
tot = 0
mons = P.monoms(); coeffs = P.coeffs()
rint=0
for mm, cfcoef in zip(mons, coeffs):
    a,b,p,q = mm  # sf^a cf^b Sh^p Ch^q
    if a%2==1: 
        continue  # odd in sf -> integrates to 0 over symmetric interval
    # integral of sin^a phi cos^(b - cf_den_power) phi ; determine den power in cf
    # den handled separately; here assume den = 2 cf^5 Sh^r ; we will divide
    tot += cfcoef * (a,b,p,q)
# Instead: do it cleanly below
print("sample kept monoms:", [mm for mm in mons if mm[0]%2==0][:5])
