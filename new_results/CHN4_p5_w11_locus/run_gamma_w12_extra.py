
import sys, time, pickle
import mpmath as mp
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M
mp.mp.dps=120
PE.F = lambda v=0, d=None: mp.mpf(v) if d is None else mp.mpf(v)/d
def coords_mpf(a, W):
    vw=PE.volreg_words_point(a,W); vwr=PE.volreg_words_point(a[::-1],W)
    T=PE.table(W); d=T.nullity; acc=[mp.mpf(0)]*d
    for src in (vw[W],vwr[W]):
        for w,pp in src.items():
            rw=T.reduce_word(w)
            for b in range(d):
                if rw[b]: acc[b]+=pp*mp.mpf(rw[b].numerator)/rw[b].denominator
    return acc
pts=pickle.load(open('gamma_extra.pkl','rb'))
res=[]
with open('gamma_w12_extra.log','w') as lg:
    for r in pts:
        a=[mp.mpf(r[0]),mp.mpf(r[1]),mp.mpf(1),mp.mpf(r[2]),mp.mpf(r[3])]
        c11=coords_mpf(a,11); c12=coords_mpf(a,12)
        CLASS=c12[9]+mp.mpf(14)/9*c12[10]
        scale=max(abs(x) for x in c12)
        m=max(abs(c12[7]),abs(c12[8]),abs(CLASS),abs(c12[11]))
        lg.write(f"a1={mp.nstr(a[0],6)} a2={mp.nstr(a[1],8)} a4={mp.nstr(a[3],8)} a5={mp.nstr(a[4],8)}\n")
        lg.write(f"  f11={mp.nstr(c11[8],4)} z3z35={mp.nstr(c11[7],4)}\n")
        lg.write(f"  z2^2z35={mp.nstr(c12[7],8)} z2z37={mp.nstr(c12[8],8)} CLASS={mp.nstr(CLASS,8)} z2118={mp.nstr(c12[11],8)}\n")
        lg.write(f"  max|W12 irred|={mp.nstr(m,8)} scale(max|c12|)={mp.nstr(scale,4)} ratio={mp.nstr(m/scale,6)}\n"); lg.flush()
        res.append({'a1':a[0],'irr':{'z2^2z35':c12[7],'z2z37':c12[8],'CLASS':CLASS,'z2118':c12[11]},'maxirr':m,'scale':scale})
pickle.dump(res, open('gamma_w12_extra.pkl','wb'))
with open('gamma_w12_extra.log','a') as lg: lg.write("DONE\n")
print("done")
