# Very quick test of Lam vs swap - are they related?
import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import mpmath as mp

mp.mp.dps = 50

def depth2(x, y):
    """Li_{2,1}(x,y) = sum_{n>m>=1} x^n y^m / (n^2 m)."""
    rho = max(abs(x), abs(x * y))
    if not 0 < rho < 1:
        return mp.mpf(0) if rho == 0 else mp.nan
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target * (1 - rho)) / mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x * y, mp.mpf(0)
    for n in range(1, nmax + 1):
        total += h / n**2
        h = x * h + x * xy_power / n
        xy_power *= x * y
    return total

r2 = mp.sqrt(2)
lam_mp = depth2(1 - r2, -r2/2)
swap_mp = depth2(1 - r2/2, -r2)

print(f"Lam  = {mp.nstr(lam_mp, 30)}")
print(f"swap = {mp.nstr(swap_mp, 30)}")

# Check the reduction expressions numerically
# Lam = -L2d*l2/4 + L2k*l2/4 - L2x1*l2/2 - L3d + L3k + 2*L3x1 + l2**2*lu/4 - l2*lu**2/4 + 3*z3/4
l2_val = mp.log(2)
lu_val = mp.log(1 + r2)
z3_val = mp.zeta(3)
L2d_val = mp.polylog(2, 2 - r2)
L2k_val = mp.polylog(2, r2 - 1)
L2x1_val = mp.polylog(2, 1 - r2)
L3d_val = mp.polylog(3, 2 - r2)
L3k_val = mp.polylog(3, r2 - 1)
L3x1_val = mp.polylog(3, 1 - r2)

lam_reduction_val = (-L2d_val*l2_val/4 + L2k_val*l2_val/4 - L2x1_val*l2_val/2 - L3d_val + L3k_val + 2*L3x1_val + l2_val**2*lu_val/4 - l2_val*lu_val**2/4 + 3*z3_val/4)

print(f"Lam reduction = {mp.nstr(lam_reduction_val, 30)}")
print(f"Diff = {mp.nstr(lam_mp - lam_reduction_val, 10)}")