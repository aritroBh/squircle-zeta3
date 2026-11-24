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

z3 = cc.h((-1,-1,0))
li3i = cc.h((-1,-1,2))
T3 = (2*arb.pi()*acb(0,1))**3
r = 24*z3 + 256*li3i + T3
print('relation residual (acb):', r.real.mid().str(5,radius=False), r.imag.mid().str(5,radius=False))
print('radii:', r.real.rad().str(3,radius=False), r.imag.rad().str(3,radius=False))
print('z3 radius:', z3.real.rad().str(3,radius=False), ' li3i radius:', li3i.real.rad().str(3,radius=False))
print('|r| ~ ', abs(r).str(5))
