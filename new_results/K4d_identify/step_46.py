import mpmath as mp, json
from mpmath import mpf
mp.mp.dps=70
WORK='research/08_theory/K4d_identify'
import sys; sys.path.insert(0,WORK)
from pslq_lib import pslq
SB=[mpf(s) for s in json.load(open(WORK+'/modelB_h96o.json'))['Sdk']]
qB=SB[1]/SB[0]
print("q_B =", mp.nstr(qB,60))
rt3=mp.sqrt(3)
# S-unit term values
tv = {
 'ln2': mp.log(2), 'ln3': mp.log(3), 'ln13': mp.log(13),
 'ln(2-r3)': mp.log(2-rt3),
 'ln((4r3-3)/13)': mp.log((4*rt3-3)/13),
 'ln((4-r3)/13)': mp.log((4-rt3)/13),
 'ln(r3-3/2)': mp.log(rt3-mpf(3)/2),
 'ln r3': mp.log(rt3),
}
for k,v in tv.items(): print("  %s = %s"%(k, mp.nstr(v,30)))
def try_pslq(tag,target,basis,tol=mpf('1e-52')):
    names=list(basis); x=[target]+[basis[n] for n in names]
    rel,resid=pslq(x,tol=tol,maxsteps=900)
    if rel is None: print("   %s: none"%tag); return None
    h=max(abs(r) for r in rel)
    print("   %s: h=%.2e resid=%.1e"%(tag,float(h),float(resid)))
    print("      rel=%s"%(rel,))
    return rel,resid
# Try: q_B vs {1} + log S-units
b1={'1':mpf(1)}; b1.update(tv)
try_pslq("{1,S-units}", qB, b1)
# smaller: maybe only a few logs matter
b2={'1':mpf(1),'ln(2-r3)':tv['ln(2-r3)'],'ln((4r3-3)/13)':tv['ln((4r3-3)/13)'],'ln((4-r3)/13)':tv['ln((4-r3)/13)'],'ln(r3-3/2)':tv['ln(r3-3/2)']}
try_pslq("{1,4tv-logs}", qB, b2)
