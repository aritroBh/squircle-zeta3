import mpmath as mp, json
from mpmath import mpf
mp.mp.dps=75
WORK='research/08_theory/K4d_identify'
import sys; sys.path.insert(0,WORK)
from pslq_lib import pslq
SB=[mpf(s) for s in json.load(open(WORK+'/modelB_h96o.json'))['Sdk']]
qB=SB[1]/SB[0]
rt3=mp.sqrt(3)
L = {
 'ln2':mp.log(2),'ln3':mp.log(3),'ln13':mp.log(13),
 'ln(2-r3)':mp.log(2-rt3),'ln(2+r3)':mp.log(2+rt3),
 'ln(4-r3)':mp.log(4-rt3),'ln(4+r3)':mp.log(4+rt3),
 'ln(5/2-r3)':mp.log(mpf(5)/2-rt3),'ln(5-2r3)':mp.log(5-2*rt3),
 'ln(r3-3/2)':mp.log(rt3-mpf(3)/2),
}
def try_pslq(tag,target,basis,tol=mpf('1e-56')):
    names=list(basis); x=[target]+[basis[n] for n in names]
    rel,resid=pslq(x,tol=tol,maxsteps=1200)
    if rel is None: print("   %s: none"%tag); return
    h=max(abs(r) for r in rel)
    star = ' <<<' if h<mpf('1e7') else ''
    print("   %s: h=%.2e resid=%.1e%s"%(tag,float(h),float(resid),star))
    if h<mpf('1e8'): print("      rel=%s"%(rel,))
def B(**kw): return kw
try_pslq("A:{1,ln2,ln3,ln(2-r3),ln(4-r3),ln(4+r3)}", qB, B(**{'1':mpf(1),'ln2':L['ln2'],'ln3':L['ln3'],'a':L['ln(2-r3)'],'b':L['ln(4-r3)'],'c':L['ln(4+r3)']}))
try_pslq("B:{1,r3,ln2,ln3,ln(2-r3),ln(4-r3)}", qB, B(**{'1':mpf(1),'r3':rt3,'ln2':L['ln2'],'ln3':L['ln3'],'a':L['ln(2-r3)'],'b':L['ln(4-r3)']}))
try_pslq("C:{1,ln(2-r3),ln(5/2-r3),ln(5-2r3)}", qB, B(**{'1':mpf(1),'a':L['ln(2-r3)'],'d':L['ln(5/2-r3)'],'e':L['ln(5-2r3)']}))
try_pslq("D:{1,ln2,ln3,ln13,ln(2-r3)}", qB, B(**{'1':mpf(1),'ln2':L['ln2'],'ln3':L['ln3'],'ln13':L['ln13'],'a':L['ln(2-r3)']}))
try_pslq("E:{1,ln2,ln3,ln(2-r3),ln(4-r3)}", qB, B(**{'1':mpf(1),'ln2':L['ln2'],'ln3':L['ln3'],'a':L['ln(2-r3)'],'b':L['ln(4-r3)']}))
try_pslq("F:{1,ln(2-r3),ln(4-r3),ln(r3-3/2)}", qB, B(**{'1':mpf(1),'a':L['ln(2-r3)'],'b':L['ln(4-r3)'],'f':L['ln(r3-3/2)']}))
