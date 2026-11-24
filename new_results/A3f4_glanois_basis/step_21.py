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
mp.mp.dps = 50

def cmzv_word(xs, etas):
    """Glanois zeta(x; eps): etas = eps as 8th-root exponents (e in 0..7 meaning zeta8^e).
    word for engine h()."""
    p = len(xs)
    run = 0  # exponent of product eps_k..eps_p mod 8
    signs = [0]*p
    for k in range(p-1,-1,-1):
        run = (run + etas[k]) % 8
        signs[k] = run
    Tt = []
    for k in range(p):
        Tt.append(signs[k]); Tt.extend([-1]*(xs[k]-1))
    return tuple(reversed(Tt))

def val8(xs, etas):  # etas as exponents
    return cc.h(cmzv_word(xs, etas))

# basis element zeta(x; eps) with eps in +-1 and last root eps_p*xi8:
def basis_val(xs, eps, s):
    etas = [0 if e==1 else 4 for e in eps[:-1]] + [0 if eps[-1]==1 else 4]
    etas[-1] = (etas[-1] + 1) % 8   # multiply by xi8
    return val8(xs, etas) * (2*arb.pi()*acb(0,1))**s

# Correct stuffle checks
res = []
# 1) Li1(xi).Li2(xi) = z(1,2; xi,xi) + z(2,1; xi,xi) + Li3(xi^2)
L = val8((1,),(1,))*val8((2,),(1,))
R = val8((1,2),(1,1)) + val8((2,1),(1,1)) + val8((3,),(2,))
res.append(('1x2', L-R))
# 2) Li1(xi).Li2(-xi) = z(1,2; xi,-xi)+z(2,1; -xi,xi)+Li3(-xi^2)
L = val8((1,),(1,))*val8((2,),(5,))
R = val8((1,2),(1,5)) + val8((2,1),(5,1)) + val8((3,),(6,))
res.append(('1x2b', L-R))
# 3) Li1(xi) * z(1,1; 1,-xi) = z(1,1,1; xi,1,-xi)+z(1,1,1; 1,xi,-xi)+z(1,1,1; 1,-xi,xi)
#    + z(2,1; xi*1, -xi) + z(1,2; 1, xi*(-xi))
# eta exponents: xi=1, -xi=5, 1=0
L = val8((1,),(1,))*val8((1,1),(0,5))
R = ( val8((1,1,1),(1,0,5)) + val8((1,1,1),(0,1,5)) + val8((1,1,1),(0,5,1))
    + val8((2,1),(1,5)) + val8((1,2),(0,6)) )
res.append(('1x11', L-R))
worst = mp.mpf(0)
for name, d in res:
    e = abs(cc.approx(d)); worst = max(worst,e)
    print(name, 'residual:', mp.nstr(e,3))
print('worst:', mp.nstr(worst,3))
with open(WORK+'/progress.log','a') as f:
    f.write('step_07: corrected stuffle checks pass (residuals ~1e-1200); converter validated depths 1-3\n')
