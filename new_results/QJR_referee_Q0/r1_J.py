# Referee QJR: independent check of Lemma J (closed form of J).
import sympy as sp, mpmath as mp
q,u,r=sp.symbols('q u r',positive=True)
R=sp.sqrt(q**2+r**2)
beta=sp.sqrt(1+q**2)
# (1) substitution u=r+R
r_u=(u**2-q**2)/(2*u); R_u=(u**2+q**2)/(2*u)
assert sp.simplify(R_u**2-(q**2+r_u**2))==0 and sp.simplify(r_u+R_u-u)==0
drdu=sp.diff(r_u,u); assert sp.simplify(drdu-(u**2+q**2)/(2*u**2))==0
orig=(1/(1+r+R)-1/(1+q))/(r*(1+r))
trans=((1/(1+u)-1/(1+q))/(r_u*(1+r_u)))*drdu
g=-2*(u**2+q**2)/((1+q)*(1+u)*(u+q)*(u**2+2*u-q**2))
print('(1) transformed integrand == g :',sp.simplify(trans-g)==0)
# (2) partial fractions with symbol B for beta, B^2=1+q^2
B=sp.symbols('B',positive=True)
gB=-2*(u**2+q**2)/((1+q)*(1+u)*(u+q)*(u+1-B)*(u+1+B))
claim=(2/((q**2-1)*(u+1))-2*q/((q**2-1)*(u+q))
       -2*(B-1)/(B*(1+q)*(B+q-1))/(u+1-B)+2*(B+1)/(B*(1+q)*(B+1-q))/(u+1+B))
d=sp.together(gB-claim); num=sp.expand(sp.numer(d)).subs(B**2,1+q**2)
num=sp.expand(num.subs(B**2,1+q**2)); num=sp.rem(sp.Poly(num,B),sp.Poly(B**2-1-q**2,B)) if num!=0 else 0
print('(2) partial fractions exact (mod B^2=1+q^2):',num==0 or sp.expand(num.as_expr())==0)
res=[2/(q**2-1),-2*q/(q**2-1),-2*(B-1)/(B*(1+q)*(B+q-1)),2*(B+1)/(B*(1+q)*(B+1-q))]
s=sp.together(sum(res)); n=sp.expand(sp.numer(s))
n=sp.rem(sp.Poly(n,B),sp.Poly(B**2-1-q**2,B)).as_expr()
print('    residues sum to 0:',sp.expand(n)==0)
# (3) residue simplifications claimed in the proof
c3=sp.together(2*(B-1)/(B*(1+q)*(B+q-1))-(1/B-1/(1+q))); n3=sp.rem(sp.Poly(sp.expand(sp.numer(c3)),B),sp.Poly(B**2-1-q**2,B)).as_expr()
c4=sp.together(2*(B+1)/(B*(1+q)*(B+1-q))-(1/B+1/(1+q))); n4=sp.rem(sp.Poly(sp.expand(sp.numer(c4)),B),sp.Poly(B**2-1-q**2,B)).as_expr()
print('(3) residue simplifications:',sp.expand(n3)==0,sp.expand(n4)==0)
print('    (q+1)^2-beta^2=2q:',sp.simplify((q+1)**2-beta**2-2*q)==0,
      ' (q+1+beta)^2=2(q+beta)(1+beta):',sp.simplify(sp.expand((q+1+beta)**2-2*(q+beta)*(1+beta)))==0,
      ' coeff log(2q):',sp.simplify(-2*q/(1-q**2)-1/(1+q)+(1+q)/(1-q**2))==0)
# (4) independent: differentiate F with true beta
F=sum(c*sp.log(u-uk) for c,uk in zip([res[0],res[1],res[2],res[3]],[-1,-q,-1+B,-1-B])).subs(B,beta)
print('(4) F\'-g == 0 (true beta):',sp.simplify(sp.diff(F,u)-g)==0)
Jcf=(2*sp.log(1+q)-(1+q)*sp.log(2*q))/(1-q**2)-(sp.asinh(q)+sp.asinh(1/q))/beta
# J=-F(q) vs closed form, numerically exact check at rationals via high-precision + symbolic positivity
mp.mp.dps=60
for qq in ['0.7071067811865475244008443621048490392848','0.72','0.75','0.85','0.95','0.999','0.3','1.7']:
    qv=mp.mpf(qq)
    f=lambda rr:(1/(1+rr+mp.sqrt(qv**2+rr**2))-1/(1+qv))/(rr*(1+rr))
    Jq=mp.quad(f,[0,qv/4,qv,1,4,40,mp.inf])
    mp.mp.dps=90
    cf=sp.N(Jcf.subs(q,sp.Rational(qq)),95)
    mF=-sp.N(F.subs(u,q).subs(q,sp.Rational(qq)),95)
    print('  q=',qq[:8],' J(quad)=',mp.nstr(Jq,25),' |quad-closed|=',mp.nstr(abs(Jq-mp.mpf(str(cf))),3),' |-F(q)-closed|=',mp.nstr(abs(mp.mpf(str(mF))-mp.mpf(str(cf))),3)); mp.mp.dps=60
# exact: -F(q) - closed form with logs as symbols (log(q+1-b)+log(q+1+b)=log 2q, ratio identity)
