t = r'''
# step_05a.py -- numerical verification of the differential formula (weight 2, fast)
import mpmath as mp
mp.mp.dps = 40

def I2(a0, a1, a2, a3):
    # int_{a0}^{a3} dt1/(t1-a1) int_{a0}^{t1} dt2/(t2-a2), a0=0, all distinct, convergent
    f = lambda t: mp.log((t - a2)/(a0 - a2)) / (t - a1)
    return mp.quad(f, [a0, a3])

def I1(a0, a1, a2):
    return mp.log((a2 - a1)/(a0 - a1))

a0, a1, a2, a3 = mp.mpf(0), mp.mpf('0.31'), mp.mpf('-0.17'), mp.mpf('0.8')
val = I2(a0, a1, a2, a3)
print("I =", mp.nstr(val, 20))

# cyclic convention: prev_1 = a3, next_2 = a0
# dI = I(a0;a2;a3) dlog((a3-a1)/(a2-a1)) + I(a0;a1;a3) dlog((a1-a2)/(a0-a2))
#      ^ i=1: [prev_1 - a_1]/[a_1 - next_1] = (a3-a1)/(a2-a1)
#      ^ i=2: [prev_2 - a_2]/[a_2 - next_2] = (a1-a2)/(a0-a2)
h = mp.mpf('1e-18')
names = ['a0','a1','a2','a3']
for k in range(4):
    pp = [a0, a1, a2, a3]; pp[k] += h
    dnum = (I2(*pp) - val) / h
    Ia = I1(a0, a2, a3)  # omit a1
    Ib = I1(a0, a1, a3)  # omit a2
    # derivative of log terms wrt param k
    def dlog_term(pp, k, num_idx, den_idx):
        # d/dp_k log((p_num - p_a)/(p_a - p_den))... handle generically below
        return 0
    # explicit: term1 = dlog((a3-a1)/(a2-a1)) * Ia ; term2 = dlog((a1-a2)/(a0-a2)) * Ib
    def d_log(rat_fun, k):
        num, den = rat_fun  # functions of the params
        # d log(num) - d log(den) wrt k-th param
        def dlogf(f):
            f0 = f(a0, a1, a2, a3)
            pp2 = [a0, a1, a2, a3]; pp2[k] += h
            return (f(*pp2) - f0) / h / f0
        return dlogf(num) - dlogf(den)
    t1 = Ia * d_log((lambda a0,a1,a2,a3: a3-a1, lambda a0,a1,a2,a3: a2-a1), k)
    t2 = Ib * d_log((lambda a0,a1,a2,a3: a1-a2, lambda a0,a1,a2,a3: a0-a2), k)
    pred = t1 + t2
    print(f"d/d{names[k]}: numeric {mp.nstr(dnum,10)}  formula {mp.nstr(pred,10)}  err {mp.nstr(abs(dnum-pred),3)}")
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_05a.py","w").write(t)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_05a.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=600)
print(r.stdout); print("ERR:", (r.stderr or "")[-1500:])
