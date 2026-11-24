
import os, json
import mpmath as mp
mp.mp.dps = 60
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs.json"))
tau = [mp.mpf(x) for x in h['eta']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = mp.arg(z0); thf = float(th)
wser = [mp.mpf(x) for x in json.load(open(wd+"/w_series_K120.json"))['w']]
def mlsq(ks, ncorr, aco, beta=1.5):
    y = [aco[k]/mp.factorial(k)*R**k for k in ks]
    nrow=2*ncorr
    S=[[mp.mpf(0)]*nrow for _ in range(nrow)]; T=[mp.mpf(0)]*nrow
    for i,k in enumerate(ks):
        bas=[]
        for j in range(ncorr):
            bas += [k**(-mp.mpf(beta)-j)*mp.cos(thf*k), k**(-mp.mpf(beta)-j)*mp.sin(thf*k)]
        for a in range(nrow):
            T[a]+=bas[a]*y[i]
            for b in range(nrow): S[a][b]+=bas[a]*bas[b]
    return mp.lu_solve(mp.matrix(S), mp.matrix(T))
print("w-series fits (beta=-3/2):")
for k0,nc in [(60,4),(80,5),(100,5)]:
    x = mlsq(list(range(k0,K+1)), nc, wser)
    Cw = mp.sqrt(x[0]**2+x[1]**2); pw = mp.atan2(x[1],x[0])
    print(f"  k0={k0} nc={nc}: C_w={mp.nstr(Cw,15)} phi_w={mp.nstr(pw,13)}")
C1 = mp.sqrt(3*mp.pi/8)
print("C_w/C_tau =", mp.nstr(mp.sqrt(x[0]**2+x[1]**2)/C1,10))
# theory: C_w = C_T * (2/sqrt3) ... C_T for T: C_tau*(sqrt3/2) (since S/S1=(2/sqrt3)(1-t)^2 T):
# B_T strength D_T; B_S/S1 = (2/sqrt3) D_T; B_w = B_T/T0 = D_T/(sqrt3/2): so C_w = C_tau*(2/sqrt3)/(sqrt3/2)?? 
# C_tau = sqrt(3 pi R/8) corresponds to B_{S/S1}. C_T = C_tau*sqrt3/2 (T vs S/S1: T = (sqrt3/2)(1-t)^-2 S/S1)
# C_w = C_T/T0 = C_tau*(sqrt3/2)/(sqrt3/2) = C_tau ???
print("C_tau*(2/sqrt3)/(sqrt3/2)*(1/1) =", mp.nstr(C1*(2/mp.sqrt(3))/(mp.sqrt(3)/2),10))
print("C_tau/2 =", mp.nstr(C1/2,12), " C_tau*2/sqrt3 =", mp.nstr(C1*2/mp.sqrt(3),12))
