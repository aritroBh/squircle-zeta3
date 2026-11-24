
import mpmath as mp
import pickle
mp.mp.dps = 80
pi=mp.pi; L2=mp.log(2)
D = pickle.load(open('/tmp/Nvals.pkl','rb'))
INT = D['INT']; N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']

# asymptotic: INT[m] ~ -2 int_0^inf (1/s) Re[(1+ i s/pi)^m] L ds
def INTasym(m):
    f = lambda s: (1/s)*mp.re((1+1j*s/pi)**m)*mp.log(1+mp.e**(-s))
    return -2*mp.quad(f, [0, mp.inf])

print("m : INT[m]                    INTasym[m]                  ratio")
for m in range(5, 43):
    a = INTasym(m)
    print(f"{m:2d} {mp.nstr(INT[m],10):>15} {mp.nstr(a,10):>15}  {mp.nstr(INT[m]/a,8)}")

# and check N_k vs sum
print()
for k in [37,38,39,40]:
    print(k, mp.nstr(N[k],8), " sum INT:", mp.nstr(sum(INT[0:k]),8))
