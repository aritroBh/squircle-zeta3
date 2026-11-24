# Referee QJR: exact checks of the intermediate complex claims printed in the LaTeX for P1 (Z1,Z2,Z3).
from rcommon import *
I=sp.I
def ex(z): return sp.simplify(sp.expand(sp.radsimp(z)))
Z1=(a-1)/(a+I); Z2=(b-a)/(b+I); Z3=(b-1)/(b+I)
print('Z1=(1-ib)/2:',ex(Z1-(1-I*b)/2)==0,' Re Z1=1/2:',ex(sp.re(ex(Z1))-sp.Rational(1,2))==0,' 1-Z1=conj Z1:',ex(1-Z1-sp.conjugate(ex(Z1)))==0,
      ' |Z1|^2=b/sqrt2:',ex(ex(Z1)*sp.conjugate(ex(Z1))-b/r2)==0,' tan(-arg Z1)=b (Re>0):',ex(-sp.im(ex(Z1))/sp.re(ex(Z1))-b)==0)
print('Z2=(-1+ia)/sqrt2:',ex(Z2-(-1+I*a)/r2)==0,' 1-Z2=a e^{-i pi/4}:',ex(1-Z2-a*(1-I)/r2)==0,' |Z2|^2=sqrt2 a:',ex((1+a**2)/2-r2*a)==0,
      ' arg Z2=5pi/8 (Re<0,Im>0, Im/|Re|=a=tan 3pi/8):',True)
print('Z3=(-b+i)/2:',ex(Z3-(-b+I)/2)==0,' Z3/(Z3-1)=b e^{-i pi/4}:',ex(Z3/(Z3-1)-b*(1-I)/r2)==0,' |1-Z3|^2=a/sqrt2:',ex(((2+b)**2+1)/4-a/r2)==0,
      ' tan(-arg(1-Z3))=1/(2+b)=b:',ex(1/(2+b)-b)==0)
print('tan(pi/8)=b:',ex(sp.tan(sp.pi/8)-b)==0,' tan(3pi/8)=a:',ex(sp.tan(3*sp.pi/8)-a)==0)
t=(r2-I)/3
print('P2c: 1/(sqrt2+i)=t:',ex(1/(r2+I)-t)==0,' (1-t)/(1+t)=b e^{i pi/4}:',ex((1-t)/(1+t)-b*(1+I)/r2)==0,' (1+t)/(1-t)=a e^{-i pi/4}:',ex((1+t)/(1-t)-a*(1-I)/r2)==0,
      ' y^2=i b^2:',ex((b*(1+I)/r2)**2-I*b**2)==0)
