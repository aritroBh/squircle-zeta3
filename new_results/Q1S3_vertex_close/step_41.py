import mpmath as mp, importlib, q3_eval
importlib.reload(q3_eval)
mp.mp.dps=60
ge = q3_eval.GEval(60, h='0.012', T='4.6', gln=5)
def panel_sum(f):
    acc=mp.mpc(0); prev=mp.mpf(0)
    for x1 in ge.xs:
        h2=(x1-prev)/2; mid=(x1+prev)/2
        s=mp.mpc(0)
        for gx,gw in ge.gl:
            s+=gw*f(mid+h2*gx)
        acc+=h2*s; prev=x1
    return acc
print("exp resid:", mp.nstr(abs(panel_sum(mp.exp)-(mp.e-1)),4))
b=mp.sqrt(2)-1
f=lambda t: mp.log1p(-t/(1j*b))/t
# integral over first half and second half separately
acc0=mp.mpc(0); prev=mp.mpf(0)
for k,x1 in enumerate(ge.xs):
    h2=(x1-prev)/2; mid=(x1+prev)/2
    s=sum(gw*f(mid+h2*gx) for gx,gw in ge.gl)
    acc0+=h2*s; prev=x1
    if k==len(ge.xs)//2: mid1=acc0
print("half1:", mp.nstr(mid1,15))
# compare exact half integral via series
v = sum((1j/b)**n*mp.mpf(1)/(n*(n+1))*(ge.xs[len(ge.xs)//2])**(n+1) for n in range(1,400))
print("half1 exact series:", mp.nstr(v,15), " resid:", mp.nstr(abs(mid1-v),4))
