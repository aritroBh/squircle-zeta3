
import mpmath as mp, json
mp.mp.dps = 60
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
tau = [mp.mpf(x) for x in h1['eta']]
d   = [mp.mpf(x) for x in h1['d']]
print(" k    log10|tau_k|   log10|d_k|    log10(k! A3^k)  log10(k! Ad^k)")
A3 = 2/mp.sqrt(mp.pi**2 + mp.log(2)**2); Ad = mp.mpf('1.28621645')
for k in [10,20,30,40,50,60,70,76,80]:
    lt = mp.log10(abs(tau[k])); ld = mp.log10(abs(d[k]))
    pt = mp.log10(mp.factorial(k)) + k*mp.log10(A3)
    pd = mp.log10(mp.factorial(k)) + k*mp.log10(Ad)
    print(f"{k:3d}  {mp.nstr(lt,8)}  {mp.nstr(ld,8)}  {mp.nstr(pt,8)}  {mp.nstr(pd,8)}")
print()
print("sign pattern d_k, k=60..80:", ''.join('+' if d[k]>0 else '-' for k in range(60,81)))
print("sign pattern tau_k, k=60..80:", ''.join('+' if tau[k]>0 else '-' for k in range(60,81)))
print()
print("d_k raw, k=74..80:")
for k in range(74,81): print(k, mp.nstr(d[k],6))
