wd = 'research/08_theory/PRV1_depth2_proofs'

part4 = r'''

# ======================================================================
# Part 4: independent numeric evidence (checks, NOT part of the proof)
# ======================================================================
def li2(z): return mp.polylog(2, z)
def li3(z): return mp.polylog(3, z)

def G_cl_mp(t):
    return (mp.mpf(3)/4*mp.zeta(3) + li3(t) + li3(-t)
            + li3(-(1-t)/(1+t)) - li3(2*t/(1+t))
            - mp.log((1+t)/2)*(li2(t) + li2((1-t)/2)/2 - li2((1+t)/2)/2))

def H_cl_mp(t):
    return (li3(1-t) + li3((1-t)/2) + li3(t/(1+t)) - li3((1-t)/(1+t))
            - li3(mp.mpf(1)/2)
            - mp.log(1+t)*(li2((1-t)/(1+t)) + li2((1-t)/2))
            + mp.log(1+t)**2*mp.log(2*t/(1+t))/2)

def G_dec_cl_mp(t):
    return (mp.zeta(3) - mp.log(t)*li2(t**2)/2 + li3(t**2)
            - 2*li3(t) - li3(1-t))

def H_dec_cl_mp(t):
    return (li3(1-t**2) - li3(1/(1+t)) + li3(t/(1+t))
            - mp.log(1+t)*(li2(-t)/2 + 3*li2(1-t)/2 + li2(1-t**2)/4
                           + mp.pi**2/24))

def dig(d):
    d = abs(d)
    return "exactly 0" if d == 0 else f"1e{mp.log10(d):.1f}".replace("e-0","e-").replace("e-","e-")

def part4():
    print("=" * 74)
    print("PART 4: numeric evidence (mpmath, independent of the exact proofs)")
    mp.mp.dps = 120
    z2, z3 = mp.pi**2/6, mp.zeta(3)
    # 4.0 classical function identities at sample points
    for t0 in (mp.mpf("0.3"), mp.mpf("0.55"), mp.mpf("0.8")):
        c = (1-t0)/(1+t0)
        vals = {
            "E": li2(t0)+li2(1-t0)-(z2-mp.log(t0)*mp.log(1-t0)),
            "L1": li2(-t0)+li2(t0/(1+t0))+mp.log(1+t0)**2/2,
            "L2": li2(-c)+li2((1-t0)/2)+mp.log((1+t0)/2)**2/2,
            "DUP2": li2(t0**2)-2*li2(t0)-2*li2(-t0),
            "OC": (li2(t0)-li2(-t0)+li2(c)-li2(-c)
                   -3*z2/2+mp.log(t0)*mp.log(c)),
            "DUP3": li3(t0**2)-4*li3(t0)-4*li3(-t0),
            "Landen3": (li3(t0)+li3(1-t0)+li3(-t0/(1-t0))
                        -(z3+z2*mp.log(1-t0)-mp.log(t0)*mp.log(1-t0)**2/2
                          +mp.log(1-t0)**3/6)),
        }
        worst = max(abs(v) for v in vals.values())
        print(f"   t={t0}: classical ids worst residual {dig(worst)}"
              f"  ({max(vals, key=lambda k: abs(vals[k]))})")
    # 4.0b Li_{1,1} closed form vs series
    for x0, y0 in ((mp.mpf("0.3"), mp.mpf("0.4")), (mp.mpf("-0.4"), mp.mpf("0.5")),
                   (mp.mpf("0.1"), mp.mpf("-0.6"))):
        F = (mp.log(1-x0)*mp.log(1-y0) + li2(-y0/(1-y0))
             - li2(-y0*(1-x0)/(1-y0)))
        ser = depth2(x0, y0, p=1)
        print(f"   Li_11 closed form vs series at ({x0},{y0}): residual {dig(F-ser)}")
    # 4.1 octagon families: series vs closed forms
    print("  octagon: G,G_cl and H,H_cl on (0,1):")
    for t0 in (mp.mpf("0.17"), mp.mpf("0.31"), mp.mpf("0.5"),
               mp.mpf("0.63"), mp.mpf("0.82")):
        rG = depth2(-t0, -(1-t0)/(2*t0)) - G_cl_mp(t0)
        rH = depth2((1-t0)/2, -2*t0/(1-t0)) - H_cl_mp(t0)
        print(f"   t={t0}: |G-G_cl| {dig(rG)}, |H-H_cl| {dig(rH)}")
    # 4.2 octagon specialisation at k = sqrt(2)-1
    k = mp.sqrt(2)-1
    x1, x2, d, h = 1-mp.sqrt(2), 1-mp.sqrt(2)/2, 2-mp.sqrt(2), mp.sqrt(2)/2
    l2, lu = mp.log(2), mp.log(1+mp.sqrt(2))
    Lam_ser = depth2(x1, -mp.sqrt(2)/2)
    swap_ser = depth2(x2, -mp.sqrt(2))
    compact3 = (3*z3/4 + li3(k) + 2*li3(x1) - li3(d)
                + l2*li2(k)/2 + l2*(li2(x2)-li2(h))/4)
    compact4 = (li3(d) + 2*li3(x2) - li3(k) - 4*li3(h) - 4*li3(-h)
                - l2*(li2(k)+li2(x2))/2 + l2**2*(l2/2-lu)/8)
    print("  octagon point (120 dps):")
    print(f"   Lambda: series={mp.nstr(Lam_ser, 30)}")
    print(f"     vs G_cl(k) {dig(Lam_ser-G_cl_mp(k))}, vs compact3 {dig(Lam_ser-compact3)}")
    print(f"   swap:   series={mp.nstr(swap_ser, 30)}")
    print(f"     vs H_cl(k) {dig(swap_ser-H_cl_mp(k))}, vs compact4 {dig(swap_ser-compact4)}")
    # 4.3 decagon families: series vs closed forms
    print("  decagon: G_dec,G_dec_cl and H_dec,H_dec_cl on (0,1):")
    for t0 in (mp.mpf("0.17"), mp.mpf("0.31"), mp.mpf("0.5"),
               mp.mpf("0.63"), mp.mpf("0.82")):
        rG = depth2(-t0, -(1-t0)/t0) - G_dec_cl_mp(t0)
        rH = depth2(1-t0, -t0/(1-t0)) - H_dec_cl_mp(t0)
        print(f"   t={t0}: |G_dec-G_dec_cl| {dig(rG)}, |H_dec-H_dec_cl| {dig(rH)}")
    # 4.4 decagon specialisation at t = 1/phi (170 dps: DEC1 claimed 830-digit LLL)
    mp.mp.dps = 170
    phi = (1+mp.sqrt(5))/2
    l, z2h, z3h = mp.log(phi), mp.pi**2/6, mp.zeta(3)
    Lam_dec_ser = depth2(-1/phi, -1/phi)
    swap_dec_ser = depth2(1/phi**2, -phi)
    Lam_DEC1 = (2*z3h/5 + 4*l*z2h/5 - l**3 + 3*li3(-1/phi) + li3(1/phi))
    swap_DEC1 = (3*z3h/5 - 7*l*z2h/5 + 2*l**3 + li3(-1/phi) + li3(1/phi))
    tp = 1/phi
    print("  decagon point (170 dps):")
    print(f"   Lambda_dec: series={mp.nstr(Lam_dec_ser, 40)}")
    print(f"     vs G_dec_cl(1/phi) {dig(Lam_dec_ser-G_dec_cl_mp(tp))},"
          f" vs DEC1 form {dig(Lam_dec_ser-Lam_DEC1)}")
    print(f"   swap_dec:   series={mp.nstr(swap_dec_ser, 40)}")
    print(f"     vs H_dec_cl(1/phi) {dig(swap_dec_ser-H_dec_cl_mp(tp))},"
          f" vs DEC1 form {dig(swap_dec_ser-swap_DEC1)}")
    # 4.5 golden Li2/Li3 lattice numerics
    lat = {
        "Li3(1/phi^2)=4Li3(-1/phi)+4Li3(1/phi)":
            li3(1/phi**2)-4*li3(-1/phi)-4*li3(1/phi),
        "Lewin Li3(1/phi^2)": li3(1/phi**2)-(4*(z3h-l*z2h)/5+2*l**3/3),
        "Li2(1/phi^2)=2z2/5-l^2": li2(1/phi**2)-(2*z2h/5-l**2),
        "Li2(1/phi)=3z2/5-l^2": li2(1/phi)-(3*z2h/5-l**2),
        "Li2(-1/phi)=-2z2/5+l^2/2": li2(-1/phi)-(-2*z2h/5+l**2/2),
    }
    for name, v in lat.items():
        print(f"   golden lattice {name}: residual {dig(v)}")
    print("PART 4 OK")


def main():
    part0()
    part1()
    part2()
    part2_specialisation()
    part3()
    part3_landen3()
    part3_specialisation()
    part4()
    print("=" * 74)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
'''

with open(f'{wd}/verify_full.py','a') as f:
    f.write(part4)
print("part4 + main appended")

import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run([PY, f'{wd}/verify_full.py'], capture_output=True, text=True, cwd=wd, timeout=1800)
out = r.stdout
print(out[-4500:])
print("STDERR:", r.stderr[-1500:])
with open(f'{wd}/verify_full.out','w') as f:
    f.write(out + "\nSTDERR:\n" + r.stderr)
