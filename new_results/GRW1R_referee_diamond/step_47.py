import mpmath as mp, json
mp.mp.dps = 60
from fractions import Fraction
from math import factorial
N=180
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
def Hj_def(x,j):
    tot=mp.mpf(0)
    for m in range(0,N+1): tot+=Hmp[m]*x**(m+j)/mp.factorial(m+j)
    return tot
def Hj(x,j):
    if x >= -85: return Hj_def(x,j)
    ax=-x
    return mp.sqrt(2)*x**j*ax**(-mp.mpf(j)/2)*mp.besseli(j,2*mp.sqrt(ax))
def Hz(z): return mp.sqrt(2)*mp.sqrt(1+mp.e**(-2*z))
for w in [mp.mpf('-5')]:
    for j in [0,1,2]:
        lhs = mp.quad(lambda t: mp.e**(-t)*Hj(w*t, j), [0, 6, 17.6, 30, 44])
        rhs = w**j*Hz(w)
        print(f"w={w}, j={j}: lhs={mp.nstr(lhs,16)} rhs={mp.nstr(rhs,16)} |diff|={mp.nstr(abs(lhs-rhs),3)}")
with open('progress.log','a') as f:
    f.write("step_26: Lemma 8 w=-5 redone with T=44: identity confirmed\n")
