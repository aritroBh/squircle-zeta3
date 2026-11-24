
import os, json, math
import mpmath as mp
mp.mp.dps = 60
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
h = json.load(open(wd+"/SRF3_hi_coeffs_K120.json"))
tau = [mp.mpf(x) for x in h['eta']]; dco = [mp.mpf(x) for x in h['d']]
K = h['K']
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0); th = float(mp.arg(z0))

# final fit table (the proper evidence): fitted model residuals
def mlsq(ks, ncorr, aco):
    y = [aco[k]/mp.factorial(k)*R**k for k in ks]
    nrow=2*ncorr
    S=[[mp.mpf(0)]*nrow for _ in range(nrow)]; T=[mp.mpf(0)]*nrow
    for i,k in enumerate(ks):
        bas=[]
        for j in range(ncorr):
            bas += [k**(-mp.mpf(1.5)-j)*mp.cos(th*k), k**(-mp.mpf(1.5)-j)*mp.sin(th*k)]
        for a in range(nrow):
            T[a]+=bas[a]*y[i]
            for b in range(nrow): S[a][b]+=bas[a]*bas[b]
    x = mp.lu_solve(mp.matrix(S), mp.matrix(T))
    return x
for nm, aco in [("tau",tau),("d",dco)]:
    ks = list(range(90,121))
    x = mlsq(ks, 5, aco)
    resid = []
    for k in ks:
        mod = sum((x[2*j]*mp.cos(th*k)+x[2*j+1]*mp.sin(th*k))*k**(-1.5-j) for j in range(5))
        resid.append(abs(aco[k]/mp.factorial(k)*R**k - mod)/abs(aco[k]/mp.factorial(k)*R**k))
    C = mp.sqrt(x[0]**2+x[1]**2); phi = mp.atan2(x[1],x[0])
    print(f"{nm}: C={mp.nstr(C,17)} phi={mp.nstr(phi,15)} max rel resid (k=90..120)={mp.nstr(max(resid),3)}")

# how many digits of C, phi confirmed
C_pred = mp.sqrt(3*mp.pi*R/8); phi_pred = th-mp.pi  # phi here = atan2 convention: -(pi-th/2)? careful
x = mlsq(list(range(90,121)), 5, tau)
C = mp.sqrt(x[0]**2+x[1]**2); phi = mp.atan2(x[1],x[0])
print("\nC err digits:", -mp.log10(abs(C-C_pred)/C_pred))
# phi convention: c_k = C cos(k th - phi_atan2) with phi_atan2 = atan2(x1,x0) = -2.2476163 = th/2 - pi
print("phi err:", abs(phi-(th/2-mp.pi)))
