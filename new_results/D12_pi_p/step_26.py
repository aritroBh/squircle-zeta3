
import mpmath as mp
import pickle
mp.mp.dps = 100
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi=mp.pi; L2=mp.log(2); z3=mp.zeta(3); z5=mp.zeta(5)

N5form = (L2**6/720 + pi**2/24*(1-L2+L2**2/2-L2**3/6) + z3/4*(1-L2+L2**2/2)
          + pi**4/2880*(1-L2) - 7*pi**4/1920 - z3*pi**2/48 + z5/8)
print("N5 closed form vs table:", mp.nstr(abs(N5form-N[5]),3))
print("N5 =", mp.nstr(N5form,40))

# leftover L6 and a_6
L6 = -7*pi**4/480 - z3*pi**2/12 + z5/2
a6form = pi**2/6 + z3 + pi**4/720 + L6
r = {k: 4*(N[k-1]-(-L2)**k/mp.factorial(k)) for k in range(3,42)}
a = {m: sum(L2**j/mp.factorial(j)*r[m-j] for j in range(0,m-2)) for m in range(3,42)}
print("a_6 form vs numeric:", mp.nstr(abs(a6form-a[6]),3), "  a_6 =", mp.nstr(a[6],30))
# equivalently a_6 = zeta(2)+zeta(3)+pi^4/720 - 7pi^4/480 - z3 pi^2/12 + z5/2
#                = zeta(2)+zeta(3)+zeta(5)/2 - z3*pi^2/12 - 13 pi^4/1440
alt = pi**2/6+z3+z5/2-z3*pi**2/12-13*pi**4/1440
print("alt form diff:", mp.nstr(abs(alt-a[6]),3))
