# step_08.py -- TASK 3 (p=1 expansion) and TASK 4 (Euclidean comparison)
import mpmath as mp, pickle
mp.mp.dps = 110
N = pickle.load(open('/tmp/Nvals.pkl','rb'))['N']; KMAX=44
f=[mp.mpf(1)]+[N[k] for k in range(KMAX+1)]
# (a) task variable eta = 1/q - 1 = -h :  pi_q/4 = 1 + sum_{k>=1} (-1)^k N_{k-1} eta^k
print("task-3 series in eta=1/q-1: coefficients (-1)^k N_{k-1}:")
for k in range(1,13): print(f"  {(-1)**k}*N_{k-1} = {mp.nstr((-1)**k*N[k-1],35)}")
# (b) alternative variable eta=q-1: h=eta/(1+eta) -> M_n
M=[mp.mpf(0)]*(KMAX+2); M[0]=mp.mpf(1)
for n in range(1,KMAX+2):
    M[n]=sum(f[j]*(-1)**(n-j)*mp.binomial(n-1,j-1) for j in range(1,n+1))
print("M_n (pi_{1+eta}/4) for n=1..15:", [mp.nstr(M[n],12) for n in range(1,16)])
print("growth |M_n|^(1/n):", [mp.nstr(abs(M[n])**(1/n),8) for n in [10,15,20,25,30,35,39]])
# verify against direct evaluation (optimal truncation)
def pi_cosh(p):
    p=mp.mpf(p); h=1/p
    g=lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(g,[0,mp.inf])
for eta in [mp.mpf('0.02'),mp.mpf('0.05'),mp.mpf('0.1')]:
    true=pi_cosh(1+eta)/4; s=mp.mpf(0); best=(0,mp.mpf(1))
    for n in range(0,42):
        s+=M[n]*eta**n; err=abs(s-true)
        if err<best[1]: best=(n,err)
    print(f"eta={eta}: best trunc n={best[0]} err={mp.nstr(best[1],3)} (pi/h with h=eta/(1+eta): {mp.nstr(mp.pi*(1+eta)/eta,6)})")
# TASK 4: Euclidean constants
print("\nEuclidean M0 = log2 - sqrt2*log(1+sqrt2) =", mp.nstr(mp.log(2)-mp.sqrt(2)*mp.log(1+mp.sqrt(2)),40))
print("our N0 = -log2 =", mp.nstr(-mp.log(2),40))
