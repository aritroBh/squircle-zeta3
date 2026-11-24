
import pickle, time
import mpmath as mp
mp.mp.dps = 110
sqrt5 = mp.sqrt(5); phi = (1+sqrt5)/2
kap1 = 1/phi; kap4 = -1/phi**2
parsed = pickle.load(open("parsed_mons_dec.pkl","rb"))
need = set()
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    for n in cu:
        need.add((n,a,b,c,d))

def integ(n,a,b,c,d):
    def f(x):
        v = x**n
        if a: v *= mp.log(x)**a
        if b: v *= mp.log(1-x)**b
        if c: v *= mp.log(1+kap1*x)**c
        if d: v *= mp.log(1+kap4*x)**d
        return v
    return mp.quad(f, [0,1], method='tanh-sinh')

t0=time.time()
ITAB = {}
for i,key in enumerate(sorted(need)):
    ITAB[key] = integ(*key)
    if (i+1)%150==0: print("  %d/%d (%.1fs)"%(i+1,len(need),time.time()-t0))
print("computed %d integrals in %.1fs"%(len(ITAB),time.time()-t0))
pickle.dump({k: str(v) for k,v in ITAB.items()}, open("itab_dec110.pkl","wb"))

# assemble P5/L = 1200 * sum (p+q sqrt5) (log2)^e (logphi)^f I(n;a,b,c,d)
L2v = mp.log(2); LFv = mp.log(phi)
acc = mp.mpf(0); acc5 = mp.mpf(0)
contrib = {'local': mp.mpf(0),'single-1': mp.mpf(0),'single-4': mp.mpf(0),'mixed': mp.mpf(0)}
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else (('single-1' if c>0 else 'single-4') if (c>0 or d>0) else 'local')
    constfac = L2v**e * LFv**f
    for n,(p,q) in cu.items():
        term = (mp.mpf(str(p)) + mp.mpf(str(q))*sqrt5)*constfac*ITAB[(n,a,b,c,d)]
        acc += term
        if sum(mon)==5:
            acc5 += term; contrib[cls] += term
gold = pickle.load(open("p5_gold.pkl","rb"))
P5L_gold = mp.mpf(gold['P5overL'])
print("\nP5/L assembled =", mp.nstr(1200*acc, 60))
print("P5/L engine    =", mp.nstr(P5L_gold, 60))
print("diff =", mp.nstr(abs(1200*acc - P5L_gold), 4))
print("\nweight-5 part W5/L =", mp.nstr(1200*acc5, 60))
for k in contrib: print("  ", k, "=", mp.nstr(1200*contrib[k], 50))
pickle.dump({'W5L_110': str(1200*acc5), 'P5L_assembled': str(1200*acc),
             'blocks': {k: str(1200*v) for k,v in contrib.items()}}, open("w5L_num.pkl","wb"))
print("saved w5L_num.pkl")
