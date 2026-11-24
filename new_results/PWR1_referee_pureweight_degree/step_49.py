import pwr1_engine as E, sympy as sp, json, time
t0=time.time()
tot3, pure3, B3 = E.total_log(3)
print("n=3 engine through tau^10 in", round(time.time()-t0,1),"s")
gam_ok = all(sp.expand(tot3.coeff(k)).coeff(E.GAM,1)==0 for k in range(11))
w_ok = all(E.weight_check(pure3.coeff(k), k) for k in range(11))
print("gamma-free:", gam_ok, "pure weight:", w_ok)
json.dump({str(k): str(pure3.coeff(k)) for k in range(11)}, open('w_exact_n3_tau10.json','w'))
for k in range(2, 9):
    print(f"w_{k}(3) = {pure3.coeff(k)}")
