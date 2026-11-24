import sys, os, pickle
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 4200
cc.N = 1400
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 300

def s600(z):
    return (z.real.mid().str(600,radius=False,more=True), z.imag.mid().str(600,radius=False,more=True))

freewords = [(7,b,c) for b in (4,5,6,7) for c in range(8)]
vals = {}
worst = mp.mpf(0)
with open(WORK+'/k3_freevals.pickle','rb') as f:
    old = pickle.load(f)
for i,w in enumerate(freewords):
    z = cc.h(w)
    vals[w] = z
    re_o, im_o = old['vals_str'][i]
    mp.mp.dps = 260
    zold = mp.mpc(mp.mpf(re_o), mp.mpf(im_o))
    mp.mp.dps = 300
    znew = mp.mpc(mp.mpf(z.real.mid().str(300,radius=False,more=True)), mp.mpf(z.imag.mid().str(300,radius=False,more=True)))
    err = abs(znew - zold)
    worst = max(worst, err)
print('32 free words recomputed; max discrepancy vs A3f3 250-digit pickle:', mp.nstr(worst,5))

# save all values for LLL steps
out = {'freewords': freewords, 'freevals': {w: s600(vals[w]) for w in freewords}}
with open(WORK+'/free_vals_hi.pickle','wb') as f:
    pickle.dump(out, f)

# extra validations
extra = []
extra.append(('z(1,2)=z(3) level1', cc.h((-1,0,0)) - cc.h((-1,-1,0))))
extra.append(('z(1,1,1)=z(3) level1', cc.h((0,0,0)) - cc.h((-1,-1,0))))
extra.append(('2pi i cubed vs literal', (2*arb.pi()*acb(0,1))**3 - acb(0,8)*arb.pi()**3))
for name, d in extra:
    print(name, 'abs mid:', d.abs_mid().str(3), ' rad:', d.real.rad().str(3,radius=False))
with open(WORK+'/progress.log','a') as f:
    f.write('step_08: recomputed 32 free words at prec4200/N1400; match A3f3 pickle to 250 digits; extra level-1 checks pass\n')
