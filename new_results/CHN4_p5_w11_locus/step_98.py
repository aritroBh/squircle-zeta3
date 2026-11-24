import subprocess
w12f = r'''
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
pts=pickle.load(open('gamma_final_pts.pkl','rb'))
res=[]
with open('gamma_final.log','w') as lg:
    for r in pts:
        a=[mp.mpf(r[0]),mp.mpf(r[1]),mp.mpf(1),mp.mpf(r[2]),mp.mpf(r[3])]
        c8=coords_mpf(a,8); c10=coords_mpf(a,10); c11=coords_mpf(a,11); c12=coords_mpf(a,12)
        CLASS=c12[9]+mp.mpf(14)/9*c12[10]
        scale=max(abs(x) for x in c12)
        m=max(abs(c12[7]),abs(c12[8]),abs(CLASS),abs(c12[11]))
        ctrl=max(abs(c8[3]),abs(c10[6]),abs(c11[7]),abs(c11[8]))
        lg.write(f"a1={mp.nstr(a[0],8)} a2={mp.nstr(a[1],8)} a4={mp.nstr(a[3],8)} a5={mp.nstr(a[4],8)}\n")
        lg.write(f"  ctrl(W8z35,W10z37,W11z3z35,f11) max={mp.nstr(ctrl,3)} ({mp.nstr(c8[3],2)},{mp.nstr(c10[6],2)},{mp.nstr(c11[7],2)},{mp.nstr(c11[8],2)})\n")
        lg.write(f"  z2^2z35={mp.nstr(c12[7],8)} z2z37={mp.nstr(c12[8],8)} CLASS={mp.nstr(CLASS,8)} z2118={mp.nstr(c12[11],8)}\n")
        lg.write(f"  max|W12 irred|={mp.nstr(m,8)} scale={mp.nstr(scale,5)} ratio={mp.nstr(m/scale,6)}\n"); lg.flush()
        res.append({'a1':a[0],'a2':a[1],'a4':a[3],'a5':a[4],'z2^2z35':c12[7],'z2z37':c12[8],'CLASS':CLASS,'z2118':c12[11],'maxirr':m,'scale':scale,'ratio':m/scale,'ctrl':ctrl})
pickle.dump(res, open('gamma_final.pkl','wb'))
with open('gamma_final.log','a') as lg: lg.write("DONE\n")
print("done")
'''
open('run_gamma_final.py','w').write(w12f)
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_gamma_final.py'],
                      stdout=open('run_gamma_final.out','w'), stderr=subprocess.STDOUT)
open('run_gamma_final.pid','w').write(str(proc.pid))
print("launched final W12 batch pid", proc.pid)
