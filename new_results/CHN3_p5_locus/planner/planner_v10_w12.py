# Planner (W4): does the weight-<=10 irreducible-free surface V10 = {C4 = Q7' = 0} of the 5-chain survive weights 11-12?
# (1) Points of V10: a = (1,1,1,a4,a5) solving C4 = Q7' = 0 (C4.txt, Q7p.txt), by Newton at 260 digits from CHN3's
#     certified isolating values (a4 ~ 1.2937, a5 ~ 0.0705) and the second positive point (a4 ~ 5.1967, a5 ~ 6.5304).
# (2) Evaluate the irreducible coordinates of c_W, W = 8..12, at those points with CHN3's validated point engine
#     (point_engine.py; exact vs forms10 at W<=10), run here in 250-digit mpf arithmetic. The m=0 (Gamma-product) part
#     carries no irreducible coordinates, so it is omitted. Controls: the W=8 zeta(3,5) and W=10 zeta(3,7) coordinates
#     must vanish at these points (they are C4 and Q7' times nonzero factors).
# Weight-12 irreducible content mod products: coef(zeta(3,9)) + (14/9) coef(zeta(5,7)) [zeta(5,7) = (14/9) zeta(3,9) mod
# products, CHN1] and coef(zeta(2,1,1,8)); also the products zeta2^2 zeta(3,5), zeta2 zeta(3,7).
import sys, os, time
import mpmath as mp
import sympy as sp
HERE = '<project>/research/08_theory/CHN3_p5_locus'
os.chdir(HERE); sys.path.insert(0, HERE)
OUT = HERE + '/planner/planner_v10_w12.out'
def log(s):
    with open(OUT, 'a') as f: f.write(time.strftime('%H:%M:%S ') + s + '\n')
import point_engine as PE, mzv_exact as M
mp.mp.dps = 260
a1, a2, a3, a4, a5 = sp.symbols('a1 a2 a3 a4 a5')
C4 = sp.sympify(open('C4.txt').read()); Q7 = sp.sympify(open('Q7p.txt').read())
c4 = sp.lambdify((a4, a5), C4.subs({a1: 1, a2: 1, a3: 1}), 'mpmath'); q7 = sp.lambdify((a4, a5), Q7.subs({a1: 1, a2: 1, a3: 1}), 'mpmath')
pts = []
for g in [(mp.mpf('1.293668784620243'), mp.mpf('0.07051999939180706')), (mp.mpf('5.1967'), mp.mpf('6.5304'))]:
    sol = mp.findroot([lambda x, y: c4(x, y), lambda x, y: q7(x, y)], g)
    x, y = sol[0], sol[1]
    log(f'V10 point: a4 = {mp.nstr(x, 30)}, a5 = {mp.nstr(y, 30)}; |C4| = {mp.nstr(abs(c4(x, y)), 3)}, |Q7p| = {mp.nstr(abs(q7(x, y)), 3)}')
    pts.append((mp.mpf(1), mp.mpf(1), mp.mpf(1), x, y))
# mpf version of the point engine: replace Fraction coercion of the point by mpf
PE.F = lambda v=0, d=None: mp.mpf(v) if d is None else mp.mpf(v) / d
mp.mp.dps = 250
LAB = {8: {3: 'z(3,5)'}, 9: {}, 10: {5: 'z2 z(3,5)', 6: 'z(3,7)'}, 11: {7: 'z3 z(3,5)', 8: 'z(3,3,5)'},
       12: {7: 'z2^2 z(3,5)', 8: 'z2 z(3,7)', 9: 'z(3,9)', 10: 'z(5,7)', 11: 'z(2,1,1,8)'}}
for W in (8, 10, 11, 12):
    log(f'basis W={W}: ' + ', '.join('*'.join(','.join(map(str, w)) for w in mon) for mon in M.BASIS[W]))
for ip, a in enumerate(pts):
    t0 = time.time()
    vw = PE.volreg_words_point(a, 12); vwr = PE.volreg_words_point(a[::-1], 12)
    log(f'point {ip}: words built ({time.time()-t0:.0f}s)')
    for W in (8, 10, 11, 12):
        T = PE.table(W); d = T.nullity
        acc = [mp.mpf(0)] * d
        for src in (vw[W], vwr[W]):
            for w, pp in src.items():
                rw = T.reduce_word(w)
                for b in range(d):
                    if rw[b]: acc[b] += pp * mp.mpf(rw[b].numerator) / rw[b].denominator
        scale = max(abs(x) for x in acc)
        irr = {LAB[W][b]: acc[b] for b in LAB[W]}
        msg = '; '.join(f'{k} = {mp.nstr(v, 8)}' for k, v in irr.items())
        if W == 12:
            cls = acc[9] + mp.mpf(14) / 9 * acc[10]
            msg += f';  CLASS z(3,9)+(14/9)z(5,7) = {mp.nstr(cls, 8)}'
        log(f'point {ip} W={W} (max |coord| {mp.nstr(scale, 4)}): {msg}  ({time.time()-t0:.0f}s)')
log('DONE')
