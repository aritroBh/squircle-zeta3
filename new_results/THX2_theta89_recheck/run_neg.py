# negative controls: the verifier must REJECT a perturbed target / perturbed certificate
from common import *
Dtarget = pickle.load(open(T+'Dtilde_relations.pkl', 'rb'))
for (W, m) in [(7, 1), (8, 2)]:
    data, B = load_rows_meta(W)
    h = norm_h(pickle.load(open(T+f'cert_D{W}_{m}.pkl', 'rb')), len(B))
    coeffs = [F(c) for c in Dtarget[f"{W},{m}"]]
    def V_of(cs):
        V = Dtilde_vec(W, m)
        for c, ns in zip(cs, ZB[W]): V = vsub(V, zprod(ns), c)
        return V
    print(W, m, "|V| terms:", len(V_of(coeffs)), " genuine:", verify_certificate(h, V_of(coeffs), data, B, False)[0])
    c2 = list(coeffs); c2[0] += F(1, 10**6)
    print("   target coeff perturbed by 1e-6 ->", verify_certificate(h, V_of(c2), data, B, False)[0])
    h2 = list(h); i = next(i for i, x in enumerate(h2) if x); h2[i] += F(1, 10**9)
    print("   certificate entry perturbed by 1e-9 ->", verify_certificate(h2, V_of(coeffs), data, B, False)[0])
