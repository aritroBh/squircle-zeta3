import mpmath as mp
mp.mp.dps = 120
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"

def M2_typed(cc):
    lam = mp.log(2); lp = mp.log(1+cc); lm = mp.log(1-cc); l0 = mp.log(cc)
    q = (1-cc**2)/(2*cc)
    Phi = mp.mpf(1)/2 - q*mp.atanh(cc)
    z2 = mp.pi**2/6
    val = (z2*Phi - mp.polylog(2, cc) - mp.polylog(2, -cc)
        + (1+cc)*mp.polylog(3, cc) + (1-cc)*mp.polylog(3, -cc)
        + ((1-cc)/cc)*mp.polylog(3, 1-cc) - ((1+cc)/cc)*mp.polylog(3, 1/(1+cc))
        - 2*cc*mp.polylog(3, (1+cc)/2)
        + (1-cc)*(mp.polylog(3, (1-cc)/(1+cc)) + mp.polylog(3, -(1-cc)/(1+cc)))
        - 2*cc*mp.polylog(3, 2*cc/(1+cc))
        + mp.polylog(2, (1+cc)/2)*((1+cc)*lp - (1-cc)*lm)
        + q*(mp.polylog(2, cc) - mp.polylog(2, -cc))*(lp+lm)
        - (1-cc)*lm**2*lp + (1-cc)*lm*lp**2 + ((4*cc**2-cc+1)/(6*cc))*lp**3
        + q*l0*lm**2 - (1-cc)*l0*lm*lp - ((1+cc**2)/(2*cc))*l0*lp**2
        + lam*((1-cc)*lm**2 - (1-cc)*lm*lp - cc*lp**2)
        + lam**2*((1-cc)/2)*(lp-lm)
        - mp.pi**2*(((1-cc)/(6*cc))*lm + ((1+cc)/(6*cc))*lp)
        + (mp.mpf(7)/4 + 2*cc)*mp.zeta(3) - (cc/6)*mp.pi**2*lam + (cc/3)*lam**3)
    return val

for tag, alpha in [("pi_2", mp.pi/2), ("2pi_3", 2*mp.pi/3), ("pi_3", mp.pi/3)]:
    st = [l for l in open(G1 + f'values_{tag}.txt') if l.startswith('M2 =')][0]
    ref = mp.mpf(st.split('=')[1].strip())
    cc = mp.cos(alpha/2)
    diff = abs(M2_typed(cc) - ref)
    print(f"{tag:>7}: |typed B2' - stored M2| = {mp.nstr(diff, 3)}  ({'100+ digits agree' if diff < mp.mpf(10)**-100 else 'FAIL'})")
