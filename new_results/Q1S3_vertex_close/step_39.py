import mpmath as mp, importlib, q3_eval
importlib.reload(q3_eval)
mp.mp.dps=60
ge = q3_eval.GEval(60, h='0.012', T='4.6', gln=5)
print("n GL nodes:", len(ge.gl), " n TS:", ge.N)
# test panel quadrature on a smooth function
b=mp.sqrt(2)-1
f=lambda t: mp.log1p(-t/(1j*b))/t
# exact integral 0->1: -Li2(-ia)
a=1+mp.sqrt(2)
ex = -mp.polylog(2,-1j*a)
# sum over panels
acc=mp.mpc(0); prev=mp.mpf(0)
for x1 in ge.xs:
    h2=(x1-prev)/2; mid=(x1+prev)/2
    s=mp.mpc(0)
    for gx,gw in ge.gl:
        s+=gw*f(mid+h2*gx)
    acc+=h2*s; prev=x1
print("panels total:", mp.nstr(acc,25))
print("exact       :", mp.nstr(ex,25))
print("resid:", mp.nstr(abs(acc-ex),4))
# check coverage: does xs span to ~1?
print("first,last xs:", mp.nstr(ge.xs[0],3), mp.nstr(ge.xs[-1],3), "sum ws:", mp.nstr(sum(ge.ws),6))
