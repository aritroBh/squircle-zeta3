import sys, os, pickle
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
import mpmath as mp
mp.mp.dps = 100

with open(WORK+'/k3_freevals.pickle','rb') as f:
    old = pickle.load(f)
i0 = old['freewords'].index((7,4,0))
re_o, im_o = old['vals_str'][i0]

# Run A: mimic A3f3: roots at 200 bits, prec 900, N=800
ctx.prec = 200
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
ctx.prec = 900
cc.N = 800
cc.series.cache_clear(); cc.half.cache_clear(); cc.g_at_one.cache_clear(); cc.h.cache_clear()
zA = cc.h((7,4,0))
mp.mp.dps = 260
zold = mp.mpc(mp.mpf(re_o), mp.mpf(im_o))
mp.mp.dps = 100
zA_mp = mp.mpc(mp.mpf(zA.real.mid().str(120,radius=False,more=True)), mp.mpf(zA.imag.mid().str(120,radius=False,more=True)))
print('A3f3-mimic (roots@200bit) vs pickle:', mp.nstr(abs(zA_mp-zold),5))

# Run B: roots at 3000 bits
ctx.prec = 3000
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
cc.N = 1200
cc.series.cache_clear(); cc.half.cache_clear(); cc.g_at_one.cache_clear(); cc.h.cache_clear()
zB = cc.h((7,4,0))
zB_mp = mp.mpc(mp.mpf(zB.real.mid().str(120,radius=False,more=True)), mp.mpf(zB.imag.mid().str(120,radius=False,more=True)))
print('correct-roots vs pickle       :', mp.nstr(abs(zB_mp-zold),5))
print('B radius:', zB.real.rad().str(3,radius=False))
print('B value 60 digits:', mp.nstr(zB_mp, 60))
