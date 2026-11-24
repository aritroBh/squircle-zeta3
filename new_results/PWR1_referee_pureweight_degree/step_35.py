import pwr1_engine as E, sympy as sp, json
zvars = {f'z{k}': E.Z.get(k, sp.Symbol(f'z{k}')) for k in range(2, 13)}
ok = True
for n, f in [(2,'routeA_n2_K12.json'), (3,'routeA_n3_K12.json')]:
    j = json.load(open(f'../../10_frontier/H1_octahedron_end/{f}'))
    eta_list = j[str(n)]['eta']
    tot, pure, B = E.total_log(n)
    S = E.ser_exp(tot)  # S_n/S_n(1) in tau
    for k in range(0, 9):
        theirs = sp.sympify(eta_list[k], locals=zvars)
        same = sp.expand(S.coeff(k) - theirs) == 0
        ok &= same
        if not same: print(f"  MISMATCH n={n} tau^{k}: mine={S.coeff(k)}\n    theirs={theirs}")
    print(f"n={n}: S/S(1) tau-series exact match vs H1 routeA through tau^8: all {all(sp.expand(S.coeff(k)-sp.sympify(eta_list[k],locals=zvars))==0 for k in range(9))}")
print("ALL:", ok)
open('progress.log','a').write("step_32-33: logB(n=2,3) exact vs H1 logform_general_n.out through eta^8; S_n/S_n(1) tau-series exact vs H1 routeA jsons (n=2,3) through tau^8\n")
