import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import mpmath as mp

mp.mp.dps = 100

def depth2(x, y):
    rho = max(abs(x), abs(x * y))
    if not 0 < rho < 1:
        raise ValueError("Invalid arguments")
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target * (1 - rho)) / mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x * y, mp.mpf(0)
    for n in range(1, nmax + 1):
        total += h / n**2
        h = x * h + x * xy_power / n
        xy_power *= x * y
    return total

r2 = mp.sqrt(2)
l2_val = mp.log(2)
lu_val = mp.log(1 + r2)
z2_val = mp.pi**2/6
z3_val = mp.zeta(3)

# Lam
lam_mp = depth2(1 - r2, -r2/2)
L2d_val = mp.polylog(2, 2 - r2)
L2k_val = mp.polylog(2, r2 - 1)
L2x1_val = mp.polylog(2, 1 - r2)
L3d_val = mp.polylog(3, 2 - r2)
L3k_val = mp.polylog(3, r2 - 1)
L3x1_val = mp.polylog(3, 1 - r2)
lam_red = -L2d_val*l2_val/4 + L2k_val*l2_val/4 - L2x1_val*l2_val/2 - L3d_val + L3k_val + 2*L3x1_val + l2_val**2*lu_val/4 - l2_val*lu_val**2/4 + 3*z3_val/4

print("LAM VERIFICATION")
print(f"Li_{{2,1}}(1-sqrt(2), -sqrt(2)/2) = {mp.nstr(lam_mp, 50)}")
print(f"Reduction expression      = {mp.nstr(lam_red, 50)}")
print(f"Difference           = {mp.nstr(abs(lam_mp - lam_red), 10)}")
print(f"Digits matched       = {mp.nstr(-mp.log10(abs(lam_mp - lam_red)), 3)}")

# swap
swap_mp = depth2(1 - r2/2, -r2)
L2mr2h_val = mp.polylog(2, -r2/2)
L2r2h_val = mp.polylog(2, r2/2)
L2x2_val = mp.polylog(2, 1 - r2/2)
L3mr2h_val = mp.polylog(3, -r2/2)
L3r2h_val = mp.polylog(3, r2/2)
L3x2_val = mp.polylog(3, 1 - r2/2)
swap_red = L2mr2h_val*l2_val/2 - L2r2h_val*l2_val/2 + L2r2h_val*lu_val/2 - L2x2_val*l2_val + L2x2_val*lu_val/2 + L3d_val - L3k_val - 4*L3mr2h_val - 4*L3r2h_val + 2*L3x2_val - l2_val**2*lu_val/2 + l2_val*lu_val**2/2 + l2_val*z2_val/2 - lu_val*z2_val/2

print("\nSWAP VERIFICATION")
print(f"Li_{{2,1}}(1-sqrt(2)/2, -sqrt(2)) = {mp.nstr(swap_mp, 50)}")
print(f"Reduction expression       = {mp.nstr(swap_red, 50)}")
print(f"Difference            = {mp.nstr(abs(swap_mp - swap_red), 10)}")
print(f"Digits matched        = {mp.nstr(-mp.log10(abs(swap_mp - swap_red)), 3)}")