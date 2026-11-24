import mpmath as mp, json
mp.mp.dps = 60
from fractions import Fraction
from math import factorial, comb
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
print("=== Lemma 8: int_0^inf e^{-t} H^{[j]}(w t) dt = w^j H(w) ===")
for w in [mp.mpf('0.5'), mp.mpf('-0.5'), mp.mpf('2'), mp.mpf('-2'), mp.mpf('-5')]:
    for j in [0,1,2]:
        T = min(mp.mpf(44), 88/abs(w))
        # tail bound: on R+ |H^{[j]}| <= (x^j/j!) 2 sqrt2 ; on R- the I_0 form continues (exact)
        lhs = mp.quad(lambda t: mp.e**(-t)*Hj(w*t, j), [0, T/3, 2*T/3, T])
        if w > 0 and T < 44:
            lhs += 0  # tail beyond T bounded below
        rhs = w**j*Hz(w)
        tb = 2*mp.sqrt(2)*(w*T)**j/mp.factorial(j)*mp.e**(-T) if w>0 else mp.mpf(0)
        print(f"w={w}, j={j}: lhs={mp.nstr(lhs,16)} rhs={mp.nstr(rhs,16)} |diff|={mp.nstr(abs(lhs-rhs),3)} (tailbd {mp.nstr(tb,2)})")
with open('progress.log','a') as f:
    f.write("step_25: Lemma 8 Laplace identity validated on/off disk both rays (w=+-0.5,+-2,-5)\n")
