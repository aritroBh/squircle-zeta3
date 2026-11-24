from common import *
import sympy as sp
t0 = time.time()
Dtarget = pickle.load(open(T+'Dtilde_relations.pkl', 'rb'))
# (a) own e-table vs stored table
stored_e = pickle.load(open(T+'e_rm_ext.pkl', 'rb'))
bad = [(m, r) for m in range(5) for r in range(0, 9) if F(stored_e.get(f"{r},{m}", '0')) != e_table(m, r) and f"{r},{m}" in stored_e]
# NOTE stored key "r,m" with first index = the (gd x)^{2*first} power? determine orientation
print("e-table check, stored key 'a,b' == [x^b](gd x)^{2a} E(x):", all(F(v) == e_table(int(k.split(',')[0]), int(k.split(',')[1])) for k, v in stored_e.items() if int(k.split(',')[0]) <= 4 and int(k.split(',')[1]) <= 11), len(stored_e))
results = {}
for W in (5, 6, 7, 8, 9):
    data, B = load_rows_meta(W)
    for m in range(W//2):
        h = pickle.load(open(T+f'cert_D{W}_{m}.pkl', 'rb'))
        h = norm_h(h, len(B))
        V = Dtilde_vec(W, m)
        coeffs = [F(c) for c in Dtarget[f"{W},{m}"]]
        for c, ns in zip(coeffs, ZB[W]): V = vsub(V, zprod(ns), c)
        ok, mism = verify_certificate(h, V, data, B, crosscheck=(W <= 8 or True))
        print(f"D~_{{{W},{m}}} certificate EXACT: {ok}   (stored-vs-regenerated relation rows differing: {mism})  t={time.time()-t0:.0f}s", flush=True)
        results[(W, m)] = ok
print("ALL 16:", all(results.values()), len(results))
