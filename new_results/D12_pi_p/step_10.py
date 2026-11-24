
import mpmath as mp
import pickle
mp.mp.dps = 110
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
mp.mp.dps = 50
print("Table of N_k (pi_p = 4(1 + sum_{k>=0} N_k / p^{k+1})):")
for k in range(0, 41):
    print(f"N_{k:2d} = {mp.nstr(N[k],40)}")
