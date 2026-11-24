
import mpmath as mp
import pickle
mp.mp.dps = 80
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi = mp.pi

print("k : ratio N_{k+1}/N_k        | normalized N_k*pi^(k-1)/(2*(k-2)!)  (odd k)")
for k in range(2, 40):
    r = N[k+1]/N[k]
    if k % 2 == 1 and k>=3:
        norm = N[k]*pi**(k-1)/(2*mp.factorial(k-2))
        print(f"k={k:2d}: {mp.nstr(r,12):>18}   norm={mp.nstr(norm,12)}  (target {(-1)**((k-1)//2)})")
    else:
        print(f"k={k:2d}: {mp.nstr(r,12):>18}")

print("\nEven k subleading: N_k*pi^(k-1)/(2*(k-2)!) (expect -> 0, rate?)")
for k in range(6,40,2):
    norm = N[k]*pi**(k-1)/(2*mp.factorial(k-2))
    print(f"k={k:2d}: norm={mp.nstr(norm,10)}   norm*k={mp.nstr(norm*k,10)}")
