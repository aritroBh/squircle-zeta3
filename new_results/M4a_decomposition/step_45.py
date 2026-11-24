import mpmath as mp
from flint import acb, arb, ctx
mp.mp.dps = 300
ctx.prec = 1000

r2 = mp.sqrt(2)
args = [(2, r2-1), (3, r2-1), (3, -(r2-1)), (3, 1/r2),
        (4, mp.mpf(1)/2), (4, r2-1), (4, -(r2-1)), (4, -(r2-1)**2),
        (4, 1/r2), (4, 1-1/r2), (4, 2-r2)]
for n, z in args:
    v_mp = mp.polylog(n, z)
    zz = acb(mp.nstr(z, 280))
    v_fl = zz.polylog(acb(n))
    # convert flint to mpf
    re = v_fl.real
    mid = re.mid()
    v_fl_mp = mp.mpf(mid.str(280, radius=False, more=True))
    d = abs(v_mp - v_fl_mp)
    rad = re.rad()
    print(f"n={n} z={mp.nstr(z,8):>12}: |mp-flint| = {mp.nstr(d,3):>10}  flint rad = {rad.str(5)}")
