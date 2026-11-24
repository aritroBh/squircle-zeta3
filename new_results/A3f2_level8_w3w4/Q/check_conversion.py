"""Independent power-series/Chen numerical controls (not a certified tail bound)."""
from cmzv import *
import mpmath as mp
from flint import acb, arb, acb_series, acb_poly, ctx

ctx.threads = 1
ctx.prec = 200
ctx.cap = 320
mp.mp.dps = 50
N = ctx.cap
roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
roots[0],roots[2],roots[4],roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)


def qa(c): return acb(str(c.p))+acb(2).sqrt()*acb(str(c.q))
def mq(c): return mp.mpf(str(c.p))+mp.sqrt(2)*mp.mpf(str(c.q))
def af(c): return acb(str(c))
def approx(z): return mp.mpc(z.real.mid().str(55,radius=False,more=True),z.imag.mid().str(55,radius=False,more=True))


@cache
def series(w,kind):
    if not w: return acb_series([1])
    assert w[-1] != -1
    tail = series(w[1:],kind)
    x = w[0]
    if x == -1:
        return acb_series([0]+[tail[j]/j for j in range(1,N)])
    pole = roots[x] if kind == 0 else (acb(1) if x == 8 else 1-roots[x])
    return (tail/acb_series([-pole,1])).integral()


@cache
def half(w,kind):
    if not w: return acb(1)
    if w[-1] == -1:
        terms = dict(sh((-1,),w[:-1])); count = terms.pop(w)
        return (half(w[:-1],kind)*acb('0.5').log()
                -sum((a*half(t,kind) for t,a in terms.items()),acb(0)))/count
    return acb_poly(series(w,kind).coeffs())(acb('0.5'))


@cache
def g_at_one(w):
    # G0(w;1)=sum G1(prefix;1/2)^rev * G0(suffix;1/2).
    out = acb(0)
    for k in range(len(w)+1):
        # kind=1 has poles 1-root; root 1 becomes zero.
        rev = tuple(-1 if x == 0 else 8 if x == -1 else x for x in reversed(w[:k]))
        out += (-1)**k*half(rev,1)*half(w[k:],0)
    return out


@cache
def h(w):
    return (-1)**sum(x != -1 for x in w)*g_at_one(tuple(-1 if x == -1 else (-x)%8 for x in w))


def evaluate(d): return sum((af(c)*h(w) for w,c in d.items()),acb(0))


def check(label,got,expected):
    assert got.real.rad() < arb('1e-35') and got.imag.rad() < arb('1e-35'), (label, got)
    error = abs(approx(got)-expected)
    log(f'{label}: abs error={mp.nstr(error,8)}')
    assert error < mp.mpf('1e-30'),(label,error)
    return error


def main():
    log(f'Numerical conversion checks: precision={ctx.prec} bits; series terms={N}; tolerance=1e-30')
    atoms = dict.fromkeys(a for c,ats in source.conjecture_terms() for a in ats)
    worst = mp.mpf(0)
    for a in atoms:
        if a[0] == 'li': expected = mp.polylog(a[1],mq(a[2]))
        elif a[0] == 'log': expected = mp.log(mq(a[1]))
        elif a[0] == 'pi2': expected = mp.pi**2
        else: expected = mp.zeta(3)
        worst = max(worst,check(repr(a),evaluate(atom(a)),expected))
    log(f'Atom checks passed: {len(atoms)}; maximum error={mp.nstr(worst,8)}')
    # Validate every original GPL pullback independently on the real path using
    # mpmath ODE Taylor series would be costly; instead compare the complete W
    # vector to the original one-dimensional elementary integral, part by part.
    lhs,rhs,diff = vectors()
    r=mp.sqrt(2); a=1+r; b=r-1
    def ell(v): return mp.log((v-b)*(v+a)/(4*v))
    def emm(v): return mp.log((a-v)*(v+b)/(4*v))
    expected = [mp.quad(lambda v: ell(v)**3*(1/(6*v)-1/(3*(v+a))),[b,1,a]),
                mp.quad(lambda v: -7*ell(v)**3/(48*v)+ell(v)**2*emm(v)/(16*v),[b,1,a])]
    for j in range(2):
        check(f'W prefactor part {j}',evaluate(lhs[j]),expected[j])
        direct = mp.mpf(0)
        for c,ats in source.conjecture_terms():
            term = mp.mpf(str((c.p,c.q)[j]))
            for aa in ats:
                if aa[0]=='li': term*=mp.polylog(aa[1],mq(aa[2]))
                elif aa[0]=='log': term*=mp.log(mq(aa[1]))
                elif aa[0]=='pi2': term*=mp.pi**2
                else: term*=mp.zeta(3)
            direct+=term
        check(f'R prefactor part {j}',evaluate(rhs[j]),direct)
        log(f'Independent elementary W-R numerical part {j}: {mp.nstr(expected[j]-direct,12)}')
    # The distribution formula as literally printed in the brief is false.
    x=mp.mpf('0.25'); n=2; d=2
    correct=d**(n-1)*(mp.polylog(n,mp.sqrt(x))+mp.polylog(n,-mp.sqrt(x)))
    log(f'Distribution control: correct Li2(x)-sum={mp.nstr(mp.polylog(n,x)-correct,8)}; brief Li2(x^2)-sum={mp.nstr(mp.polylog(n,x*x)-correct,12)}')
    log('Numerical conversion controls PASS (uncertified truncation).')

if __name__ == '__main__': main()


@cache
def real_series(poles):
    if not poles: return acb_series([1])
    assert poles[-1] != Q(0)
    tail=real_series(poles[1:])
    if poles[0] == Q(0): return acb_series([0]+[tail[j]/j for j in range(1,N)])
    return (tail/acb_series([-qa(poles[0]),1])).integral()


@cache
def real_g(w,c,e):
    if not w: return acb(1)
    if c==e: return acb(0)
    def taylor(w,c,e):
        if w[-1]==c: return None
        nz=[abs((p-c)/(e-c)) for p in w if p!=c]
        if nz and min(nz)<=Q(F(3,2)): return None
        poles=tuple((p-c)/(e-c) for p in w)
        return acb_poly(real_series(poles).coeffs())(acb(1))
    value=taylor(w,c,e)
    if value is not None: return value
    value=taylor(tuple(reversed(w)),e,c)
    if value is not None: return (-1)**len(w)*value
    mid=(c+e)/2
    return sum((real_g(w[:k],mid,e)*real_g(w[k:],c,mid) for k in range(len(w)+1)),acb(0))


def check_each_gpl():
    worst=mp.mpf(0); count=0
    for (w,k),coefficient in source.w4_g_expression().items():
        for e in (A,B):
            pulled=forms_path([dlog(p) for p in w],POINTS[Q(1)],POINTS[e])
            got=evaluate(pulled); expected=real_g(w,Q(1),e)
            err=abs(approx(got-expected)); worst=max(worst,err); count+=1
            assert err<mp.mpf('1e-30'),(w,e,err)
    log(f'ALL individual original GPL endpoint conversions checked={count}; maximum error={mp.nstr(worst,10)}')
