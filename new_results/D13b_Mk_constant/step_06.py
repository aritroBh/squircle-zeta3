
import mpmath as mp
mp.mp.dps = 80

def B(z):
    return 1 - mp.exp(z)/(1+mp.sqrt(1+mp.exp(2*z)))

# branch points z_s = s*i*pi/2 ; cut where 1+e^{2z} negative real => zeta positive real
for s in [+1,-1]:
    zs = s*1j*mp.pi/2
    print("s =", s)
    for r in [mp.mpf('1e-3'), mp.mpf('1e-4')]:
        for eps in [mp.mpf('1e-12'), -mp.mpf('1e-12')]:
            z = zs + r + 1j*eps
            val = B(z)
            lhs = (val - (1 - s*1j))/mp.sqrt(r)   # divide by zeta^{1/2}, zeta = r (positive real)
            print(f"  r={mp.nstr(r,3)} eps={mp.nstr(eps,3)}: (B-(1-si))/sqrt(zeta) = {mp.nstr(lhs,12)}")
    # discontinuity across the cut
    r = mp.mpf('1e-4')
    disc = (B(zs+r+1j*mp.mpf('1e-12')) - B(zs+r-1j*mp.mpf('1e-12')))/(2*mp.sqrt(r))
    print("  Disc/(2 sqrt(r)) =", mp.nstr(disc, 12), "   (expect +/- s*sqrt(2) =", mp.nstr(s*mp.sqrt(2),6), ")")
