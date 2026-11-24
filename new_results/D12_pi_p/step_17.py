
import mpmath as mp
import pickle
mp.mp.dps = 100
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi = mp.pi; L2 = mp.log(2); z3 = mp.zeta(3)

# r_k = [h^k] R : N_{k} = (-log2)^(k+1)/(k+1)! + r_{k+1}/4
r = {}
for k in range(3, 9):
    r[k] = 4*(N[k-1] - (-L2)**k/mp.factorial(k))
print("r_3 =", mp.nstr(r[3],30), " pi^2/6 =", mp.nstr(pi**2/6,30))
print("r_4 =", mp.nstr(r[4],30), " pi^2/6+z3-pi^2log2/6 =", mp.nstr(pi**2/6+z3-pi**2*L2/6,30))
print("r_5 =", mp.nstr(r[5],30))

# a_m = [h^m](2^h R) = sum_j log2^j/j! r_{m-j}
for m in range(3,8):
    a = sum(L2**j/mp.factorial(j)*r[m-j] for j in range(0,m-2))
    print(f"a_{m} = {mp.nstr(a,30)}")
print()
print("candidate a_5 = pi^2/6+z3+pi^4/720 =", mp.nstr(pi**2/6+z3+pi**4/720,30))
