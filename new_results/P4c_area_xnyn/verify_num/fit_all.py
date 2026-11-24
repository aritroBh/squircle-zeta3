import mpmath as mp, json, time
from fit import load, basis, fit
from claims import claims, implied, rho_pred
DPS = 170
mp.mp.dps = DPS
C = claims(DPS); I = implied(DPS)
out = {}
names = [((0,0),'a0 (sqrt2)'),((1,2),'R: L^2 (-1/2)'),((1,1),'beta'),((1,0),'gamma'),((2,3),'R3: L^3 (1/3)'),((2,2),'b3'),((2,1),'c3'),((2,0),'d3'),
         ((3,4),'R4: L^4'),((3,3),'R4: L^3'),((3,2),'R4: L^2'),((3,1),'R4: L^1'),((3,0),'R4: L^0'),((4,5),'R5: L^5 (1/30?)')]
truth = {(0,0):C['a0'],(1,0):C['R'][0],(1,1):C['R'][1],(1,2):C['R'][2],(2,0):C['R3'][0],(2,1):C['R3'][1],(2,2):C['R3'][2],(2,3):C['R3'][3],(4,5):mp.mpf(1)/30}
for j in range(5): truth[(3,j)] = I['R4'][j]
data = load((1, 1))
req = [1000, 3000, 10000, 30000, 100000, 1000000]
cfgs = [('A', 10**3, 10**21, 10), ('B', 10**4, 10**21, 10), ('C', 10**5, 10**21, 9), ('D', 10**6, 10**21, 8), ('E', 10**3, 10**12, 9),
        ('F', 10**3, 10**6+1, 4)]
res = {}
for tag, a, b, K in cfgs:
    t = time.time()
    o, mx, m, p = fit(data, basis(K), a, b)
    res[tag] = dict(fit={str(k): mp.nstr(v, 60) for k, v in o.items()}, maxres=mp.nstr(mx, 3), rows=m, params=p, nmin=a, nmax=b, K=K)
    print(tag, a, b, K, m, p, mp.nstr(mx, 3), round(time.time()-t, 1), flush=True)
    for kj, nm in names:
        if kj in o: print('   ', nm, mp.nstr(o[kj], 30), 'diff', mp.nstr(o[kj] - truth[kj], 3))
# six-point fit at the requested n only: nE = sqrt2 + (c12 L^2 + beta L + gamma)/n + (c23 L^3+b3 L^2+c3 L+d3)/n^2 ... too many unknowns;
# use: a0 fixed sqrt2, R and R3 free except leading-log coefficients fixed, plus one effective R4 term c*L^4/n^3 -> 6 unknowns, 6 points
six = [(n, E) for n, E in data if n in req]
assert len(six) == 6
B6 = [(0,0),(1,2),(1,1),(1,0),(2,3),(2,2),(2,1),(2,0),(3,4)]
o6, mx6, m6, p6 = fit(six, B6, 0, 10**7, fixed={(0,0): mp.sqrt(2), (1,2): mp.mpf(-1)/2, (2,3): mp.mpf(1)/3})
print('six-point', m6, p6, mp.nstr(mx6, 3))
for kj, nm in names:
    if kj in o6: print('   ', nm, mp.nstr(o6[kj], 20), 'diff', mp.nstr(o6[kj] - truth[kj], 3))
res['six'] = dict(fit={str(k): mp.nstr(v, 40) for k, v in o6.items()}, maxres=mp.nstr(mx6, 3))
# loose-degree and half-power tests
o, mx, m, p = fit(data, basis(9, extra_deg=1), 10**3, 10**21)
print('loose deg+1', m, p, mp.nstr(mx, 3)); res['loose'] = {}
for kj in [(1,3),(2,4),(3,5)]:
    print('    extra coeff', kj, mp.nstr(o[kj], 5)); res['loose'][str(kj)] = mp.nstr(o[kj], 5)
for kj, nm in names[:8]:
    print('   ', nm, 'diff', mp.nstr(o[kj] - truth[kj], 3)); res['loose'][nm] = mp.nstr(o[kj] - truth[kj], 3)
o, mx, m, p = fit(data, basis(10, half=True), 10**3, 10**21)
print('half-power', m, p, mp.nstr(mx, 3), 'coef n^{-1/2} in nE:', mp.nstr(o[('h',0)], 5), 'd3 diff', mp.nstr(o[(2,0)] - truth[(2,0)], 3))
res['half'] = dict(coef=mp.nstr(o[('h',0)], 5), d3diff=mp.nstr(o[(2,0)] - truth[(2,0)], 3), maxres=mp.nstr(mx, 3))
json.dump(res, open('fit_rho1.json', 'w'), indent=1)
