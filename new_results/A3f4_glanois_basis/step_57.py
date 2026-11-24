import sys, os
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

li = cc.h((-1,-1,2))
im = li.imag
pim = arb.pi()**3/32
print('li.imag:', im.mid().str(110, radius=False))
print('pi^3/32:', pim.mid().str(110, radius=False))
d = im - pim
print('diff mid:', d.mid().str(10, radius=False), 'rad:', d.rad().str(3, radius=False))
