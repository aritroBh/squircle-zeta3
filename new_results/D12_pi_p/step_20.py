
import mpmath as mp
import pickle
mp.mp.dps = 110
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
# coefficient list of f(h)=pi_{1/h}/4: f[0]=1, f[n]=N_{n-1}
KMAX=44
f = [mp.mpf(1)]+[N[k] for k in range(KMAX+1)]
# h = eta/(1+eta) = sum_{t>=1} (-1)^{t+1} eta^t ; compose
M = [mp.mpf(0)]*(KMAX+2)
for n in range(0,KMAX+2):
    # M_n = sum_{j=1..n} f[j] * (-1)^{n-j} binom(n-1, j-1)
    tot = mp.mpf(0)
    for j in range(1,n+1):
        tot += f[j]*(-1)**(n-j)*mp.binomial(n-1,j-1)
    M[n] = tot if n>0 else mp.mpf(1)
mp.mp.dps = 45
print("Expansion at p=1:  pi_{1+eta} = 4 * sum M_n eta^n")
for n in range(0,16):
    print(f"M_{n:2d} = {mp.nstr(M[n],35)}")
pickle.dump({'M':M}, open('/tmp/Mvals.pkl','wb'))
