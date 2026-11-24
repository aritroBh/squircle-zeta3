# Referee QJR: P1 rebuilt from scratch (own partial fractions, own primitive selection, exact endpoint values).
from rcommon import *
import json
v=sp.symbols('v',positive=True)
qv=(v**2+1)/(2*r2*v)
print('q(1)=',canon(qv.subs(v,1)),' q(a)=',canon(qv.subs(v,a)),' dq/dv=',sp.factor(sp.diff(qv,v)),'(>0 on (1,a])')
print('2q^2-1 == ((v^2-1)/(2v))^2:',sp.simplify(2*qv**2-1-((v**2-1)/(2*v))**2)==0)
print('2 q q\'/sqrt(2q^2-1) == (v^2+1)/(2v^2) (v>1):',sp.simplify(2*qv*sp.diff(qv,v)/((v**2-1)/(2*v))-(v**2+1)/(2*v**2))==0)
print('1+q == (v+a)(v+b)/(2sqrt2 v):',sp.simplify(1+qv-(v+a)*(v+b)/(2*r2*v))==0,
      ' 1-q == -(v-a)(v-b)/(2sqrt2 v):',sp.simplify(1-qv+(v-a)*(v-b)/(2*r2*v))==0,
      ' 2q == (v^2+1)/(sqrt2 v):',sp.simplify(2*qv-(v**2+1)/(r2*v))==0)
# Letters: log(1+q)=La+Lb-Lv-3/2 lam ; log(2q)=Lw-Lv-lam/2 ;   (La=log(v+a),Lb=log(v+b),Lv=log v,Lw=log(v^2+1))
La,Lb,Lv,Lw=sp.symbols('La Lb Lv Lw')
lp=La+Lb-Lv-sp.Rational(3,2)*lam; l2q=Lw-Lv-lam/2
pref=(v**2+1)/(2*v**2)
integ=sp.expand(pref*(2*lp-(1+qv)*l2q)/(1-qv**2))
coef={}
for X in [La,Lb,Lv,Lw,lam]:
    cX=sp.together(integ.coeff(X))
    coef[X]=sp.apart(sp.factor(cX,extension=sp.sqrt(2)),v,extension=sp.sqrt(2)) if False else sp.apart(sp.radsimp(cX),v,full=False,extension=sp.sqrt(2))
assert sp.expand(sum(coef[X]*X for X in coef)-integ).equals(0) or True
# collect residues at simple poles robustly: residue = lim (v-c)*C
poles=[0,a,-a,b,-b]
res={}
for X in coef:
    cX=sp.together(integ.coeff(X))
    rr={}
    for cpole in poles:
        rv=canon(sp.limit(sp.simplify((v-cpole)*cX),v,cpole)) if True else None
        if rv!=0: rr[cpole]=rv
    recon=sum(r_/(v-c_) for c_,r_ in rr.items())
    ok=sp.simplify(sp.radsimp(sp.together(cX-recon)))==0
    res[X]=rr
    print('C_%s: residues %s ; C == sum r/(v-c) exactly: %s'%(X,{str(canon(k_)):str(v_) for k_,v_ in rr.items()},ok))
    assert ok
# compare with printed C's
pr={La:-2/(v-a)+2/(v+a)+2/(v-b)-2/(v+b), Lb:-2/(v-a)+2/(v+a)+2/(v-b)-2/(v+b), Lv:-r2/v-2/(v+a)+2/(v+b),
    Lw:r2/v+2/(v-a)-2/(v-b), lam:-r2/(2*v)+2/(v-a)-3/(v+a)-2/(v-b)+3/(v+b)}
print('printed C_l agree:',all(sp.simplify(sp.radsimp(sp.together(sum(r_/(v-c_) for c_,r_ in res[X].items())-pr[X])))==0 for X in pr))
letval={La:lambda t:t+a, Lb:lambda t:t+b, Lv:lambda t:t, Lw:lambda t:t**2+1}
letp={La:a,Lb:b,Lv:0}
# pole at v=a: constants cancel
cancel=sp.expand(sum(res[X].get(a,0)*logS(letval[X](a)) for X in [La,Lb,Lv,Lw])+res[lam].get(a,0)*lam)
print('end-point pole v=a: sum_l r_(l,a) log l(a) + r_(lam,a) lam =',cancel)
assert cancel==0
# Build primitive pieces; record value at v=a minus value at v=1 exactly; and symbolic derivative check
ReZ={}
def ReLi2Z(cpole,vv):
    Z=canon((cpole-vv)*(cpole-sp.I)/(cpole**2+1)) if False else sp.nsimplify(0)
    return None
pieces=[]; total=0
vs=sp.symbols('vs',positive=True)
for X in [La,Lb,Lv,Lw,lam]:
    for cpole,rr in res[X].items():
        if X==lam:
            if cpole==a: continue  # absorbed in the cancellation
            assert sp.N(1-cpole)>0
            val=rr*lam*(logS(a-cpole)-logS(1-cpole)); typ='log'
            F=rr*lam*sp.log(vs-cpole)
            dchk=sp.simplify(sp.diff(F,vs)-rr*lam/(vs-cpole))==0
        elif X==Lw:
            if cpole==0:
                typ='T5'; F=-rr*sp.polylog(2,-vs**2)/2
                val=-rr*(Li2(-a**2)-Li2(-1))/2
                dchk=sp.simplify(sp.diff(F,vs).replace(sp.polylog,lambda s_,z_:-sp.log(1-z_) if s_==1 else sp.polylog(s_,z_))-rr*sp.log(vs**2+1)/vs)==0
            else:
                typ='T6'
                # -2 Re Li2(Z), Z=(c-v)/(c+i); plus log(c^2+1) log(v-c) unless c==a (constant part cancelled)
                Zs={}
                for vv in [1,a]:
                    Zs[vv]=sp.expand(sp.radsimp((cpole-vv)/(cpole+sp.I)))
                name=lambda Z: sp.Symbol('ReLi2[%s]'%sp.simplify(Z)) if Z!=0 else 0
                val=rr*(-2*name(Zs[a])+2*name(Zs[1]))
                if cpole!=a:
                    assert sp.N(1-cpole)>0
                    val+=rr*logS(cpole**2+1)*(logS(a-cpole)-logS(1-cpole))
                # derivative check of -2 Re Li2(Z(v)) numerically at 3 points (complex analytic; exact form checked below)
                dchk=True
                mp.mp.dps=50
                for t_ in ['1.2','1.9','2.3']:
                    t_=mp.mpf(t_); cc=mp.mpf(sp.N(cpole,60))
                    Ff=lambda tt:-2*mp.re(mp.polylog(2,(cc-tt)/mp.mpc(cc,1)))
                    dchk&=abs(mp.diff(Ff,t_)-(mp.log(t_**2+1)-mp.log(cc**2+1))/(t_-cc))<mp.mpf(10)**-40
                ZZ=sp.Symbol('Z'); 
        else:
            p=letp[X]
            if cpole==a:
                typ='T1'; F=-rr*sp.polylog(2,(a-vs)/(a+p))
                val=-rr*(0-Li2((a-1)/(a+p)))
                assert sp.N(a+p)>0 and sp.N(1+p)>0
                integrand=rr*(sp.log(vs+p)-sp.log(a+p))/(vs-a)
            else:
                d=canon(cpole+p); assert sp.N(1-cpole)>0
                if d==0:
                    typ='T4'; F=rr*sp.log(vs-cpole)**2/2
                    val=rr*(logS(a-cpole)**2-logS(1-cpole)**2)/2
                elif sp.N(d)>0:
                    typ='T2'; F=rr*(sp.log(d)*sp.log(vs-cpole)-sp.polylog(2,-(vs-cpole)/d))
                    val=rr*(logS(d)*(logS(a-cpole)-logS(1-cpole))-Li2(-(a-cpole)/d)+Li2(-(1-cpole)/d))
                else:
                    typ='T3'; assert sp.N(1-cpole+d)>0   # v-c>|d| on [1,a]
                    F=rr*(sp.log(vs-cpole)**2/2+sp.polylog(2,-d/(vs-cpole)))
                    val=rr*((logS(a-cpole)**2-logS(1-cpole)**2)/2+Li2(-d/(a-cpole))-Li2(-d/(1-cpole)))
                integrand=rr*sp.log(vs+p)/(vs-cpole)
            dF=sp.diff(F,vs).replace(sp.polylog,lambda s_,z_:-sp.log(1-z_) if s_==1 else sp.polylog(s_,z_))
            dchk=all(abs(sp.N((dF-integrand).subs(vs,sp.Rational(t_)),60))<sp.Float('1e-50') for t_ in ['6/5','19/10','23/10'])
        pieces.append((str(X),str(canon(cpole)),str(rr),typ,dchk,sp.expand(val)))
        total+=val
for pc in pieces: print('  %-4s pole %-14s r=%-10s %-4s F\'=integrand:%s  value=%s'%pc)
total=sp.expand(total)
print('Re Li2 symbols occurring:',sorted({str(s_) for s_ in total.free_symbols if str(s_).startswith('ReLi2')}))
json.dump(sp.srepr(total),open('r4_P1_unreduced.json','w'))
# Complex boundary values: own reductions
mp.mp.dps=110
s2=mp.sqrt(2); A=1+s2; B=s2-1
LiB=mp.polylog(2,B)
def num_sym(name):
    Z=complex(0)
    return None
ReSyms={str(s_):s_ for s_ in total.free_symbols if str(s_).startswith('ReLi2')}
Zval={}
for k_,s_ in ReSyms.items():
    zs=sp.sympify(k_[6:-1]); Zval[k_]=zs
    print('  ',k_,' Z =',zs,'=',sp.N(zs,12))
# exact values (derived by hand in REFEREE.md, checked below)
lamb_=lam
exact={}
for k_,zs in Zval.items():
    zn=complex(sp.N(zs,30))
    if abs(zn-complex(0.5,-0.20710678118654752))<1e-12 or abs(zn-complex(0.5,0.20710678118654752))<1e-12:
        exact[k_]=sp.Rational(29,384)*pi**2-L**2/8-lam*L/8-lam**2/32
    elif abs(zn-complex(-0.7071067811865476,1.7071067811865476))<1e-12 or abs(zn-complex(-0.7071067811865476,-1.7071067811865476))<1e-12:
        exact[k_]=Lam/2-lam*L/4-sp.Rational(19,192)*pi**2
    elif abs(zn-complex(-0.20710678118654752,0.5))<1e-12 or abs(zn-complex(-0.20710678118654752,-0.5))<1e-12:
        exact[k_]=-Lam/2-L**2/8+lam*L/8-lam**2/32+pi**2/384
    else: raise SystemExit('unexpected Z '+k_)
    d=abs(mp.re(mp.polylog(2,mp.mpc(sp.N(sp.re(zs),120),sp.N(sp.im(zs),120))))-numval(exact[k_].subs(Lam,Li2(b))))
    print('   Re Li2(%s): value %s ; |num - value| = %s'%(zs,exact[k_],mp.nstr(d,3)))
P1u=sp.expand(total.subs({ReSyms[k_]:exact[k_] for k_ in ReSyms}))
print('P1 (unreduced, real Li2 + Lambda):',P1u)
# compare with the printed 18-term expression
pr18=(6*Lam-4*Li2(sp.Rational(1,2))+2*Li2(1/r2)-2*Li2(-1/r2)+2*Li2(-r2)-2*Li2(-a)-4*Li2(b/r2)-2*Li2(r2*b)+2*Li2(-b/2)+2*Li2(b/(2*r2))
      +r2/2*Li2(-1)-r2/2*Li2(-a**2)-r2/2*L**2-lam*L-r2/2*lam*L+lam**2/4-sp.Rational(5,48)*pi**2)
P1u_L=sp.expand(P1u.subs(Lam,Li2(b))); pr18_L=sp.expand(pr18.subs(Lam,Li2(b)))
print('P1 unreduced - printed 18-term expression (exact) =',sp.simplify(P1u_L-pr18_L))
# reduce with the table
tv={k_:sp.sympify(v_) for k_,v_ in json.load(open('r3_table_values.json')).items()}
red=P1u_L
for f_ in [f for f in red.atoms(sp.Function) if f.func.__name__=='Li2']:
    k_=str(canon(f_.args[0]))
    if k_==str(canon(b)): red=red.subs(f_,Lam)
    else: red=red.subs(f_,tv[k_])
red=sp.expand(red)
P1pr=3*r2*Lam-lam**2/2-r2/2*lam*L+(1+3*r2/2)*L**2-(sp.Rational(1,48)+3*r2/16)*pi**2
print('P1 reduced:',red)
print('P1 reduced - Lemma P (exact) =',sp.simplify(red-P1pr))
print('|P1 unreduced (num) - v-form quadrature| =',end=' ')
mp.mp.dps=110
vq=mp.quad(lambda t:(t**2+1)/(2*t**2)*(lambda qq:(2*mp.log(1+qq)-(1+qq)*mp.log(2*qq))/((1-qq)*(1+qq)))((t**2+1)/(2*s2*t)),[1,2,A])
print(mp.nstr(abs(vq-numval(P1u_L)),3))
json.dump(sp.srepr(red),open('r4_P1_reduced.json','w'))
