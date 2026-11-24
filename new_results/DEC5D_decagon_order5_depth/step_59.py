
import pickle, time
from fractions import Fraction as F
import dec5d_lib as D
from qsqrt2 import QS2

OCT5D = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed_oct = pickle.load(open(OCT5D+"parsed_mons.pkl","rb"))
kit2 = D.make_kit(QS2, D.logvec2, ('log2','log(1+sqrt2)'))
RW2 = D.make_rewrite(QS2, QS2(2,1), QS2(-1,-1))
ZERO2 = QS2(0)
t0 = time.time()
need = {}
for mon, cu in parsed_oct:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    for n in cu: need[(n,a,b,c,d)] = mon
rewrite2 = {k: RW2['monomial_integral'](*k) for k in sorted(need)}
print("octagon rewrite: %d integrals in %.1fs" % (len(rewrite2), time.time()-t0))
obs2, _, _ = D.make_obstruction(kit2)
gsymcache = {}
def Gsym2(w):
    if w not in gsymcache: gsymcache[w] = kit2['g_symbol'](kit2['ZERO'], w, kit2['ONE'])
    return gsymcache[w]
def top_symbol2(a,b,c,d,n):
    k = a+b+c+d
    out = {}
    for w, gam in rewrite2[(n,a,b,c,d)].items():
        if len(w) != k: continue
        out = kit2['add'](out, kit2['scale'](Gsym2(w), gam))
    return out
def qs2f(p,q): return QS2(F(p.numerator,p.denominator), F(q.numerator,q.denominator))
sub_sym = {'local': {}, 'single-p': {}, 'single-q': {}, 'mixed': {}}
t0 = time.time()
for mon, cu in parsed_oct:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else (('single-p' if c>0 else 'single-q') if (c>0 or d>0) else 'local')
    for n,(p,q) in cu.items():
        s = top_symbol2(a,b,c,d,n)
        if e: s = kit2['shuffle'](s, kit2['const_pow'](0, e))
        if f: s = kit2['shuffle'](s, kit2['const_pow'](1, f))
        sub_sym[cls] = kit2['add'](sub_sym[cls], kit2['scale'](s, qs2f(p,q)))
sub_sym = {k: kit2['scale'](v, QS2(960)) for k,v in sub_sym.items()}
print("octagon symbols built in %.1fs" % (time.time()-t0))
exp = {'local': (QS2(0), QS2(0)),
       'single-p': (QS2(F(-53955,8), F(4770)), QS2(F(-10215,2), F(14445,4))),
       'single-q': (QS2(F(-53955,8), F(4770)), QS2(F(-10215,2), F(14445,4))),
       'mixed':   (QS2(F(48015,4), F(-8505)), QS2(F(11115), F(-15525,2)))}
ok = True
for k in ('local','single-p','single-q','mixed'):
    o = obs2(sub_sym[k])
    good = o[0][0]==exp[k][0] and o[1][0]==exp[k][1]
    ok &= good
    print(f"  {k:9s}: ({o[0][0]}, {o[1][0]})   {'PASS' if good else 'FAIL exp '+str(exp[k])}")
tot = kit2['add'](*sub_sym.values())
oT = obs2(tot)
print("TOTAL octagon W5/u:", (oT[0][0], oT[1][0]), " expected (-1485+1035 sqrt2, 900-540 sqrt2)")
print("TOTAL", "PASS" if oT[0][0]==QS2(-1485,1035) and oT[1][0]==QS2(900,-540) else "FAIL")
pickle.dump({'sub_sym': sub_sym}, open("octagon_control_syms.pkl","wb"))
