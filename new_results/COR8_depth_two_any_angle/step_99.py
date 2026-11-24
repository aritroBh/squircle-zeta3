
import pickle
import mpmath as mp
mp.mp.dps = 80
DS = {eval(k): mp.mpf(v) for k,v in pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/series_D_v3.pkl','rb'))['DS'].items()}
L2 = mp.log(2)
basis = {'z4': mp.zeta(4), 'z3l2': mp.zeta(3)*L2, 'z2l2^2': mp.zeta(2)*L2**2, 'l2^4': L2**4,
         'Li4(1/2)': mp.polylog(4, mp.mpf('0.5')), 'z3': mp.zeta(3), 'z2l2': mp.zeta(2)*L2,
         'l2^3': L2**3, 'Li3(1/2)': mp.polylog(3, mp.mpf('0.5')), 'z2': mp.zeta(2),
         'l2^2': L2**2, 'Li2(1/2)': mp.polylog(2, mp.mpf('0.5')), 'l2': L2, '1': mp.mpf(1)}
blist = list(basis.items())
for key in [(-1,0), (0,0), (1,0), (-1,2), (0,1), (1,1), (1,2), (1,3), (2,0), (3,0)]:
    v = DS[key]
    x = [v] + [b for _, b in blist]
    rel = mp.pslq(x, tol=mp.mpf('1e-70'), maxcoeff=10**10)
    if rel:
        terms = [(rel[0], 'D'), *[(rel[i+1], blist[i][0]) for i in range(len(blist))]]
        lhs = [f"({c_})*{nm}" for c_, nm in terms[1:] if c_]
        print(f"D[{key}] = {mp.nstr(v, 20)}  relation: {rel}")
    else:
        print(f"D[{key}] = {mp.nstr(v, 20)}  NO relation")
