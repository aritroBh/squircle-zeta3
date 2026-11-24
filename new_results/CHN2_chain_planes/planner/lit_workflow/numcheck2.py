import sys; sys.dont_write_bytecode=True
from mpmath import mp, mpf, gamma, hyp3f2
mp.dps=35
def F1(s12,s23,s34,s45,s24):   # SS (3.11)
    return gamma(1+s12)*gamma(1+s23)*gamma(1+s34)*gamma(1+s45)/(gamma(1+s12+s23)*gamma(1+s34+s45))*hyp3f2(s12,1+s45,-s24,1+s12+s23,1+s34+s45,1)
def F2(s12,s23,s34,s45,s24,s13):
    return s13*s24*gamma(1+s12)*gamma(1+s23)*gamma(1+s34)*gamma(1+s45)/(gamma(2+s12+s23)*gamma(2+s34+s45))*hyp3f2(1+s12,1+s45,1-s24,2+s12+s23,2+s34+s45,1)
def dict_(u1,u2,u3,u4):
    return dict(s12=u3,s23=mpf(0),s34=u2,s45=u1+u3,s24=u4,s13=u1)
for u in [(mpf('0.3'),mpf('0.5'),mpf('0.7'),mpf('0.9')),(mpf('1.1'),mpf('0.2'),mpf('0.45'),mpf('2.3'))]:
    d=dict_(*u)
    f1=F1(d['s12'],d['s23'],d['s34'],d['s45'],d['s24']); f2=F2(d['s12'],d['s23'],d['s34'],d['s45'],d['s24'],d['s13'])
    lm=gamma(1+u[0]+u[2])*gamma(1+u[1]+u[3])/gamma(1+sum(u))
    print('u=',[mp.nstr(x,4) for x in u])
    print('   F1=%s F2=%s'%(mp.nstr(f1,20),mp.nstr(f2,20)))
    print('   F1-F2 - Gamma(1+u1+u3)Gamma(1+u2+u4)/Gamma(1+sum u) = %s'%mp.nstr(f1-f2-lm,5))
# Pi1/Pi2 check lambda_+ = F1+F2
G=lambda x: gamma(1+x)
x,y=mpf('0.37'),mpf('0.81')
for lab,u,lp in (('Pi1',(x,x,y,y),G(x)**2*G(y)**2/(G(2*x)*G(2*y))),('Pi2',(x,y,x,y),G(x)**2*G(y)**2/G(x+y)**2)):
    d=dict_(*u)
    f1=F1(d['s12'],d['s23'],d['s34'],d['s45'],d['s24']); f2=F2(d['s12'],d['s23'],d['s34'],d['s45'],d['s24'],d['s13'])
    print(lab,'F1+F2 - closed lambda_+ =',mp.nstr(f1+f2-lp,5))
