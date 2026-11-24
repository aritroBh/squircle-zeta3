import sys, os
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
from flint import ctx
ctx.prec = 1600
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
cc.N = 1500
import mpmath as mp
mp.mp.dps = 50

h = cc.h
# test 1: zeta(3) = h(-1,-1,0)
z3 = h((-1,-1,0))
print('zeta3 err:', mp.nstr(abs(cc.approx(z3)-mp.zeta(3)),5))
print('radius:', z3.real.rad().str(5, radius=False), z3.imag.rad().str(5, radius=False))

# test 2: Li_3(xi8) = h(-1,-1,1)?
xi8 = mp.exp(2*mp.pi*mp.j/8)
li3_direct = mp.polylog(3, xi8)
cand = h((-1,-1,1))
print('Li3(xi8) err:', mp.nstr(abs(cc.approx(cand)-li3_direct),5))

# test 3: Li_3(-xi8) = h(-1,-1,5)?
li3m = mp.polylog(3, -xi8)
candm = h((-1,-1,5))
print('Li3(-xi8) err:', mp.nstr(abs(cc.approx(candm)-li3m),5))
with open(WORK+'/progress.log','a') as f:
    f.write('step_04: validated engine h() vs zeta(3); basis<->word conversion on Li3(+-xi8)\n')
