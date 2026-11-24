from flint import arb
rho = arb.pi()/2
Estar = arb(5)/24
sqrt2 = arb(2).sqrt()
def sgn_xa(n): return [-1,1,1,-1][n%4]
def sgn_xb(n): return [1,1,-1,-1][n%4]
def sgn_law(n): return [-1,1,1,-1][n%4]
def beta(x,y): return x.gamma()*y.gamma()/(x+y).gamma()
def A_n(n):
    return sgn_xa(n)*sqrt2*(2/arb.pi())*rho**(arb('0.5')-n)*beta(arb('1.5'), arb(n)-arb('0.5'))
def B_n(n):
    return sgn_xb(n)*(sqrt2/2)*(2/arb.pi())*rho**(arb('1.5')-n)*beta(arb('2.5'), arb(n)-arb('1.5'))
def law2(n):
    return sgn_law(n)*rho**(-n)*arb(n)**arb('-1.5')*(1 + (3-3*arb.pi()*((-1)**n))/(8*n))
def Ebnd(n):
    return (4/arb.pi())*Estar*rho**(arb('2.5')-n)*beta(arb('3.5'), arb(n)-arb('2.5'))
def tails(n):
    rhot = arb('2.9')*rho
    return (4/arb.pi())*rhot**(arb('0.5')-n)/(arb(n)-arb('0.5')) + (2/arb.pi())*rhot**(arb('1.5')-n)/(arb(n)-arb('1.5'))
mx=None; argmx=None
for n in range(10,501):
    u = rho**(-n)*arb(n)**arb('-3.5')
    val = (abs(A_n(n)+B_n(n)-law2(n)) + Ebnd(n) + tails(n))/u
    if mx is None or val.upper() > mx: mx, argmx = val.upper(), n
print("certified sup_{10<=n<=500} (|A+B-law2|+Ebnd+tails)/(rho^-n n^-7/2) =", mx, " at n =", argmx)
# also the Ebound piece alone and |A+B-law2| piece alone
mxE=None; mxAB=None
for n in range(10,501):
    u = rho**(-n)*arb(n)**arb('-3.5')
    vE = (Ebnd(n)/u).upper(); vAB = (abs(A_n(n)+B_n(n)-law2(n))/u).upper()
    if mxE is None or vE>mxE: mxE=vE
    if mxAB is None or vAB>mxAB: mxAB=vAB
print("sup Ebnd/u =", mxE, " sup |A+B-law2|/u =", mxAB)
with open('progress.log','a') as f:
    f.write(f"step_09b: certified sup R_bound/u over 10..500 = {mx} at n={argmx}; supEb={mxE}, supAB={mxAB}\n")
