
import mpmath as mp
import pickle
mp.mp.dps = 100
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi = mp.pi; L2 = mp.log(2); z3 = mp.zeta(3)

# N4 closed form candidate:
N4form = -L2**5/120 + pi**2/24*(1 - L2 + L2**2/2) + z3/4*(1 - L2) + pi**4/2880
print("N4 form vs table:", mp.nstr(abs(N4form - N[4]),3))

# A(h): pi_{1/h}/4 = 2^{-h} (1 + A(h)), A = sum a_m h^m
KMAX=44
mp.mp.dps = 110
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
fser = [mp.mpf(1)] + [N[k] for k in range(KMAX+1)]  # coefficients of pi/4 in h: f[0]=1, f[n]=N_{n-1}
fser = fser[:KMAX+2]
# multiply by 2^h = sum (log2)^j/j! h^j
Aser = [mp.mpf(0)]*(KMAX+2)
for n in range(KMAX+2):
    Aser[n] = sum(fser[n-j] * L2**j/mp.factorial(j) for j in range(n+1))
print("A(h) = pi/4 * 2^h - coefficients a_m/4? :")
for m in range(1,12):
    print(f"a_{m} = 4*{mp.nstr(Aser[m],30)}")
