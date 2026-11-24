
import mpmath as mp
import pickle
mp.mp.dps = 80
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi = mp.pi; L2 = mp.log(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
Li4 = mp.polylog(4, mp.mpf('0.5')); Li5 = mp.polylog(5, mp.mpf('0.5'))

def pslq_report(target, basis, names, tag, maxcoeff=10**6):
    v = [target]+basis
    rel = mp.pslq(v, tol=mp.mpf('1e-70'), maxcoeff=maxcoeff, maxsteps=10000)
    if rel is None:
        print(f"{tag}: PSLQ found NOTHING")
        return None
    # rel[0]*target + sum rel[i+1]*basis[i] = 0
    resid = abs(rel[0]*target + sum(rel[i+1]*basis[i] for i in range(len(basis))))
    print(f"{tag}: relation found, residual={mp.nstr(resid,3)}")
    expr = ""
    for i in range(len(basis)):
        if rel[i+1]:
            expr += f" {'+' if rel[i+1]>0 else '-'} {abs(rel[i+1])/abs(rel[0])}*{names[i]}"
    print(f"   {tag} = -1/{rel[0]}*(...): {tag} ≈ {expr}")
    # verify
    val = -sum(rel[i+1]*basis[i] for i in range(len(basis)))/rel[0]
    print(f"   check: |{tag} - value| = {mp.nstr(abs(target-val),3)}")
    return rel

# N_4, weights <=4
basis4 = [mp.mpf(1), L2, L2**2, L2**3, L2**4, pi**2, pi**2*L2, pi**2*L2**2, pi**4, z3, z3*L2, Li4]
names4 = ['1','log2','log2^2','log2^3','log2^4','pi^2','pi^2*log2','pi^2*log2^2','pi^4','zeta3','zeta3*log2','Li4(1/2)']
r4 = pslq_report(N[4], basis4, names4, "N_4")

# N_5, weights <=5
basis5 = basis4 + [L2**5, pi**4*L2, pi**2*L2**3, z3*L2**2, z3*pi**2, Li4*L2, z5, Li5]
names5 = names4 + ['log2^5','pi^4*log2','pi^2*log2^3','zeta3*log2^2','zeta3*pi^2','Li4*log2','zeta5','Li5(1/2)']
r5 = pslq_report(N[5], basis5, names5, "N_5")
