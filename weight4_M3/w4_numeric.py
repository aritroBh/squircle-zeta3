"""NUMERICAL evidence only: fresh quadrature and the paper's stored A/B data."""
from w4_symbol import *
import mpmath as mp


def number(x):
    x = Q.c(x)
    return mp.mpf(x.p.numerator)/x.p.denominator + mp.sqrt(2)*mp.mpf(x.q.numerator)/x.q.denominator


def atom_value(a):
    if a[0] == 'log':
        return mp.log(abs(number(a[1])))
    if a[0] == 'li':
        return mp.polylog(a[1],number(a[2]))
    if a == ('pi2',):
        return mp.pi**2
    assert a == ('z3',)
    return mp.zeta(3)


def rhs_parts():
    out = [mp.mpf(0),mp.mpf(0)]
    for c,atoms in conjecture_terms():
        val = mp.fprod(atom_value(a) for a in atoms)
        out[0] += mp.mpf(c.p.numerator)/c.p.denominator*val
        out[1] += mp.mpf(c.q.numerator)/c.q.denominator*val
    return out


def direct_parts(dps):
    with mp.workdps(dps+30):
        a,b = number(A),number(B)
        def integrals(u):
            # v=b+(a-b)*u. log1p avoids cancellation in log(1-t).
            v = b+(a-b)*u
            t = (v-b)*(v+a)/(4*v)
            mt = (a-v)*(v+b)/(4*v)
            lt = mp.log(t) if t < mp.mpf('0.5') else mp.log1p(-mt)
            lm = mp.log(mt) if mt < mp.mpf('0.5') else mp.log1p(-t)
            return (lt**3/v,lt**3/(v+a),lt**2*lm/v)
        vals = [mp.quad(lambda u: integrals(u)[i]*(a-b),[0,mp.mpf('0.1'),mp.mpf('0.5'),1]) for i in range(3)]
        return [vals[0]/6-vals[1]/3,-7*vals[0]/48+vals[2]/16]


def tail_value():
    r2,b = mp.sqrt(2),number(B)
    lam,L = mp.log(2),mp.log(number(A))
    w3p = ((7*r2/32-mp.mpf(3)/16)*mp.zeta(3)+(1+r2)*mp.polylog(3,b)
           +(1-r2)*mp.polylog(3,-b)-2*r2*mp.polylog(3,1/r2)
           +(2*L-r2*lam)*mp.polylog(2,b)+(r2/12-mp.mpf(1)/8)*mp.pi**2*L
           -(1+r2)/24*mp.pi**2*lam+(1+r2)/24*lam**3+r2/4*lam**2*L
           -(1+r2)/2*lam*L**2+(5-r2)/6*L**3)
    return w3p+(mp.pi**2+6*lam**2)/24


def main():
    dps = int(sys.argv[1]) if len(sys.argv)>1 else 160
    mp.mp.dps = dps+30
    exactish = rhs_parts()
    direct = direct_parts(dps)
    print('NUMERICAL, fresh direct quadrature, requested digits:',dps,flush=True)
    print('W4 =',mp.nstr(direct[0]+mp.sqrt(2)*direct[1],65),flush=True)
    for i in range(2):
        err = direct[i]-exactish[i]
        print('coefficient part',i,'quadrature minus RHS =',mp.nstr(err,8),flush=True)
        print('constant c_'+str(i),'= residual/pi^4 =',mp.nstr(err/mp.pi**4,8),flush=True)
        assert abs(err)<mp.mpf(10)**(-dps+5)
    with mp.workdps(1250):
        p,qv = rhs_parts()
        val = p+mp.sqrt(2)*qv
        for method in 'AB':
            stored = mp.mpf((PAPER/'outputs'/f'M3_method{method}_1200.txt').read_text())
            w4 = stored-tail_value()
            err = w4-val
            print('NUMERICAL, stored method',method,'W4 minus RHS =',mp.nstr(err,8),flush=True)
            print('residual divided by pi^4 =',mp.nstr(err/mp.pi**4,8),flush=True)
            assert abs(err)<mp.mpf('1e-1195')


if __name__ == '__main__':
    main()
