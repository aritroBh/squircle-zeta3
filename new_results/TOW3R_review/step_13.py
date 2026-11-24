import sys, os, time, pickle
from fractions import Fraction as Fr
from math import comb
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work+'/tow3_copy'); os.chdir(work+'/tow3_copy')
import sympy as sp
import tow3lib as T3
d, m = T3.d, T3.m; i = sp.I
pS = sp.Symbol('p')

H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gdiv(x,y):
    den = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den, (x[1]*y[0]-x[0]*y[1])/den)
def gsub(x,y): return (x[0]-y[0], x[1]-y[1])
def qsi_div(x, y):
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    den = (gmul(C,C)[0]-2*gmul(D,D)[0], gmul(C,C)[1]-2*gmul(D,D)[1])
    numA = (gmul(A,C)[0]-2*gmul(B,D)[0], gmul(A,C)[1]-2*gmul(B,D)[1])
    numB = (gmul(B,C)[0]-gmul(A,D)[0], gmul(B,C)[1]-gmul(A,D)[1])
    return (gdiv(numA, den), gdiv(numB, den))
def qsi_gmul(x, c):
    A,B = x; return (gmul(A,c), gmul(B,c))

rB0d = 2*sp.I*(4*pS-1)/(pS+1)
def ehat(delta, s, pv):
    tot = 0
    for dpp in range(0, delta+1):
        tot += T3.tau1(delta-dpp) * T3.e(dpp, s, 2*pv+1, pv+1) * rB0d.subs(pS,pv) / (4*pv-dpp-1)
    return sp.expand((pv+1)*tot)

t0=time.time()
for delta in [5,6]:
    smax = delta//2
    bad=[]; tested=0
    for key in H:
        if key=='_meta': continue
        p,q=key
        if p!=q: continue
        dd=2*p
        if delta not in H[key] or 2*dd-4 < delta: continue
        e00 = H[key][0][dd]
        vals=[]; okk=True
        for kk in range(dd+1, dd+2+smax):   # k = 2p+1 .. 2p+1+smax
            if kk not in H[key][delta]: okk=False; break
            vals.append(qsi_div(H[key][delta][kk], e00))
        if not okk: continue
        # v_t = sum_s binom(t, s) E_s with t = k-(2p+1) = 0..smax
        E=[]
        for tt in range(smax+1):
            acc = vals[tt]
            for s in range(tt):
                term = qsi_gmul(E[s], (Fr(comb(tt,s)),Fr(0)))
                acc = (gsub(acc[0],term[0]), gsub(acc[1],term[1]))
            E.append(acc)
        for s in range(smax+1):
            want = ehat(delta, s, p)
            wr, wi = sp.Rational(sp.re(want)), sp.Rational(sp.im(want))
            got = E[s]; tested += 1
            if got[1]!=(Fr(0),Fr(0)) or got[0]!=(wr,wi):
                bad.append((key,s,got,(wr,wi)))
    print(f'depth {delta}: diagonal ACTUAL ehat vs formula: tested {tested}, bad {len(bad)} {bad[:3]}  ({time.time()-t0:.0f}s)', flush=True)

# empirical E3: Delta_pp(delta) from ACTUAL data (exceptional k=2p value minus stable extrapolation),
# feed it via direct + fold (TOW1 formula) at depth delta; check 0 for delta=5,6, p=3..10
def qsi_add(x,y): return ((x[0][0]+y[0][0],x[0][1]+y[0][1]),(x[1][0]+y[1][0],x[1][1]+y[1][1]))
def qsi_conj(x): return ((x[0][0],-x[0][1]),(x[1][0],-x[1][1]))  # conjugate of A+B s2: conj A + conj B s2
def qsi_negs(x): return ((-x[0][0],-x[0][1]),(-x[1][0],-x[1][1]))
def qsi_gmulq(x, c):  # QSI * QSI not needed; QSI * Gaussian rational
    A,B=x; return (gmul(A,c), gmul(B,c))
def tau1v(n):
    aq=sp.sqrt(2)+1; bq=sp.sqrt(2)-1
    s=n+1
    v = (-1)**s * ( (sp.I-bq)**(-(s+1)) + (sp.I+aq)**(-(s+1)) - sp.I**(-(s+1)) )
    v = sp.simplify(v)
    return v  # element of Q(i, sqrt2); represent as pair? use sympy throughout below
def tau2v(n):
    aq=sp.sqrt(2)+1; bq=sp.sqrt(2)-1
    s=n+1
    v = (-1)**s * ( (sp.I-aq)**(-(s+1)) + (sp.I+bq)**(-(s+1)) - sp.I**(-(s+1)) )
    return sp.simplify(v)
# represent QSI as sympy a+b*sqrt2 with a,b in Q(i)
def to_sym(x): return (x[0][0]+x[0][1]*sp.I) + (x[1][0]+x[1][1]*sp.I)*sp.sqrt(2)
for delta in [5,6]:
    smax=delta//2
    for pv in range(3, 11):
        dd=2*pv
        if (pv,pv) not in H or delta not in H[(pv,pv)]: continue
        # actual ehat from data (as above)
        e00 = H[(pv,pv)][0][dd]
        E={}
        vals=[]
        for kk in range(dd+1, dd+2+smax):
            vals.append(qsi_div(H[(pv,pv)][delta][kk], e00) if kk in H[(pv,pv)][delta] else None)
        if any(vv is None for vv in vals): continue
        for tt in range(smax+1):
            acc = vals[tt]
            for s in range(tt):
                term = qsi_gmulq(E[s], (Fr(comb(tt,s)),Fr(0)))
                acc = (gsub(acc[0],term[0]), gsub(acc[1],term[1]))
            E[tt]=acc
        # Delta in units of e00: rho(delta)/e00 - sum_s (-1)^s E_s  (rho(delta) = actual k=2p value)
        rhoN = qsi_div(H[(pv,pv)][delta][dd], e00)
        Dlt = rhoN
        for s in range(smax+1):
            term = qsi_gmulq(E[s], (Fr((-1)**s),Fr(0)))
            Dlt = (gsub(Dlt[0],term[0]), gsub(Dlt[1],term[1]))
        # feed: direct + fold, evaluated as sympy Q(i,sqrt2) expression (in units of e00_child = real*sqrt2;
        # we only check it vanishes, so units don't matter)
        Dsym = to_sym(Dlt)      # Delta (in units of e00_pp)
        Dc = -sp.conjugate(Dsym)  # conjugate in units of e00_pp (anti-real e00: actual conj = -conj(xi))
        feed = sp.Integer(0)
        for dp in range(0, delta+1):
            # need Delta at ALL depths dp<=delta: recompute recursively for dp
            pass
        # (compute Delta for all depths dp<=delta first)
        Ds={}; Dcs={}
        for dp in range(0, delta+1):
            smaxp = dp//2
            vals=[]
            ok=True
            for kk in range(dd+1, dd+2+smaxp):
                if kk not in H[(pv,pv)][dp]: ok=False; break
                vals.append(qsi_div(H[(pv,pv)][dp][kk], e00))
            if not ok: break
            E2={}
            for tt in range(smaxp+1):
                acc = vals[tt]
                for s in range(tt):
                    term = qsi_gmulq(E2[s], (Fr(comb(tt,s)),Fr(0)))
                    acc = (gsub(acc[0],term[0]), gsub(acc[1],term[1]))
                E2[tt]=acc
            rhoN = qsi_div(H[(pv,pv)][dp][dd], e00)
            Dlt = rhoN
            for s in range(smaxp+1):
                term = qsi_gmulq(E2[s], (Fr((-1)**s),Fr(0)))
                Dlt = (gsub(Dlt[0],term[0]), gsub(Dlt[1],term[1]))
            Ds[dp] = to_sym(Dlt)
            Dcs[dp] = -sp.conjugate(Ds[dp])
        else:
            direct = sum( pv*tau2v(delta-dp)*Ds[dp]/(4*pv-3-dp) for dp in range(0, delta+1) )
            fold = sp.Integer(0)
            for t in range(0, delta+1):
                mu = 4*pv-4-t
                g = sum( -Dcs[dp]*sp.conjugate(tau1v(t-dp))/(4*pv-3-dp) for dp in range(0, t+1) )
                binom = sp.binomial(mu-2, delta-t)
                fold += -pv * g * (-sp.I)**(-4-t) * binom * sp.I**(-6-delta)
            tot = sp.simplify(direct + fold)
            print(f'E3 empirical delta={delta} p={pv}: direct+fold (from ACTUAL data) = {tot}  (==0: {tot==0})', flush=True)
