
import pickle, time
from fractions import Fraction as F
import dec5d_lib as D

QS5 = D.QS5; ZERO = QS5(0)
kit5 = D.make_kit(QS5, D.logvec5, ('log2','logsqrt5','logphi'))
RW5 = D.make_rewrite(QS5, -D.PHI, D.PHI**2)
rewrite = pickle.load(open("rewrite_dec_exact.pkl","rb"))
parsed = pickle.load(open("parsed_mons_dec.pkl","rb"))

gsymcache = {}
def Gsym(w):
    if w not in gsymcache: gsymcache[w] = kit5['g_symbol'](kit5['ZERO'], w, kit5['ONE'])
    return gsymcache[w]
def top_symbol(a,b,c,d,n):
    k = a+b+c+d
    out = {}
    for w, gam in rewrite[(n,a,b,c,d)].items():
        if len(w) != k: continue
        out = kit5['add'](out, kit5['scale'](Gsym(w), gam))
    return out
def qs5f(p,q): return QS5(F(p.numerator,p.denominator), F(q.numerator,q.denominator))
sub_sym = {'local': {}, 'single-1': {}, 'single-4': {}, 'mixed': {}}
t0 = time.time()
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else (('single-1' if c>0 else 'single-4') if (c>0 or d>0) else 'local')
    for n,(p,q) in cu.items():
        s = top_symbol(a,b,c,d,n)
        if e: s = kit5['shuffle'](s, kit5['const_pow'](0, e))   # log2 -> slot 0
        if f: s = kit5['shuffle'](s, kit5['const_pow'](2, f))   # logphi -> slot 2
        sub_sym[cls] = kit5['add'](sub_sym[cls], kit5['scale'](s, qs5f(p,q)))
sub_sym = {k: kit5['scale'](v, QS5(1200)) for k,v in sub_sym.items()}
print("decagon weight-5 block symbols built in %.1fs" % (time.time()-t0))
for k, s in sub_sym.items():
    print("  %-9s: %d words" % (k, len(s)))
tot = kit5['add'](*sub_sym.values())
print("  total    : %d words" % len(tot))
slots_seen = set()
for s in list(sub_sym.values())+[tot]:
    for w in s:
        slots_seen.update(w)
print("slots appearing in any symbol:", sorted(slots_seen))
pickle.dump({'sub_sym': sub_sym, 'tot': tot}, open("dec_block_symbols.pkl","wb"))
