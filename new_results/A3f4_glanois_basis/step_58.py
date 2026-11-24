import sys, os
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 4200
ctx.cap = 6000
cc.N = 1400
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)

li = cc.h((-1,-1,2))
pim = arb.pi()**3/32
d = li.imag - pim
print('Im Li3(i) - pi^3/32 mid:', d.mid().str(4, radius=False), 'rad:', d.rad().str(3, radius=False))

# relation 24 zeta3 + 256 Li3(i) + T^3
z3 = cc.h((-1,-1,0))
T3 = (2*arb.pi()*acb(0,1))**3
r = 24*z3 + 256*li + T3
print('24z3+256Li3i+T3 mid:', r.real.mid().str(4,radius=False), r.imag.mid().str(4,radius=False))
print('radius:', r.real.rad().str(3,radius=False))
with open(WORK+'/progress.log','a') as f:
    f.write('step_20: ROOT CAUSE of A3f3 failure: ctx.cap=320 in check_conversion capped arithmetic at 320 bits; raised cap; Li3(i) now matches pi^3/32 and 24z3+256Li3(i)+T^3=0 holds\n')
