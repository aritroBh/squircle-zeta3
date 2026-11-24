
import os, json
import mpmath as mp
mp.mp.dps = 60
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs.json"))
tau = [mp.mpf(x) for x in h['eta']]; dco = [mp.mpf(x) for x in h['d']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = mp.arg(z0)
# w-series: log[(2-p)S/(4sqrt3)] = sum w_k tau^k: compute from tau coeffs
# (2-p)*S/S(1) = (1-2t)/(1-t) * S/S(1)
Sp = [mp.mpf(0)]*(K+1)
for k in range(K+1):
    # multiply (1-2t)/(1-t) = (1-2t)*sum t^j
    Sp[k] = tau[k] - 2*(tau[k-1] if k>=1 else 0) + sum((tau[j]-2*(tau[j-1] if j>=1 else 0)) for j in range(0,k))
# careful: (1-2t)/(1-t) * A = B => b_k = a_k - a_{k-1} - sum_{j<k} a_j*(...)... do it properly
# (1-2t)*A = (1-t)*B => b_k - b_{k-1} = a_k - 2 a_{k-1}
A = tau
B = [mp.mpf(0)]*(K+1)
B[0] = A[0]
for k in range(1,K+1):
    B[k] = A[k]-2*A[k-1]+B[k-1]
wser = [mp.mpf(0)]*(K+1)
# log B: w = log B (B[0]=1)
L = [mp.mpf(0)]*(K+1)
for k in range(1,K+1):
    s = mp.mpf(0)
    for j in range(1,k+1):
        s += j*B[j]*L[k-j]
    L[k] = (B[k]*k - s)/k/1  # (log B)' B = B' => k B_k = sum_j j L_j B_{k-j}... careful
# redo: log B = L; L' B = B'; k L_k = k B_k - sum_{j=1}^{k-1} j L_j B_{k-j}
L = [mp.mpf(0)]*(K+1)
for k in range(1,K+1):
    L[k] = B[k] - sum(j*L[j]*B[k-j] for j in range(1,k))/k
wser = L
print("w_2..w_6:", [mp.nstr(wser[k],12) for k in range(2,7)])
print("paper: -pi^2/9=", mp.nstr(-mp.pi**2/9,12), " -16z3/9=", mp.nstr(-16*mp.zeta(3)/9,12), " -37pi^4/810=", mp.nstr(-37*mp.pi**4/810,12))

# growth law evidence table: law/actual with correction series for the three series
C1 = mp.sqrt(3*mp.pi/8); C2 = mp.sqrt(3*mp.pi/4)
def law_tau(k): return -C1*mp.factorial(k)*k**(-mp.mpf('1.5'))*mp.re(z0**(mp.mpf('0.5')-k))
def law_d(k):   return +C2*mp.factorial(k)*k**(-mp.mpf('1.5'))*mp.im(z0**(mp.mpf('0.5')-k))
def law_w(k):   return None
print("\nk : tau: a_k/law_tau (with 1/k corr) : d: a_k/law_d")
for k in [60,80,100,110,120]:
    rt = tau[k]/law_tau(k); rd = dco[k]/law_d(k)
    print(f"  k={k}: {mp.nstr(rt,10)}   {mp.nstr(rd,10)}")
# w-series law: B_w has same z0 with strength D_T/T(0)... w_k/k! ~ -(2/sqrt3)*sqrt(3pi/8)*... derive: B_{log}~B_T/T(0):
# B_T ~ (3 pi i/4) e^{z0} z0^{1/2}(1-z/z0)^{1/2}: w law = B_T/T0 strength with T0=sqrt3/2:
# a_w/k! ~ -(1/sqrt pi) k^-3/2 2 Re[D_T/T0 z0^-k]... compute directly:
T0 = mp.sqrt(3)/2
DT = (3*mp.pi*1j/4)*mp.exp(z0)*z0**mp.mpf('0.5')
def law_w(k):
    return -(1/mp.sqrt(mp.pi))*mp.factorial(k)*k**(-mp.mpf('1.5'))*2*mp.re(DT/T0*z0**(-k))
print("\nw-series: w_k/law_w:")
for k in [60,80,100,120]:
    print(f"  k={k}: {mp.nstr(wser[k]/law_w(k),10)}")
json.dump({'w':[mp.nstr(wser[k],50) for k in range(K+1)]}, open(wd+"/w_series_K120.json","w"))
print("saved w_series_K120.json")
