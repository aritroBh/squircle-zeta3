import mpmath as mp
mp.mp.dps = 100
import numpy as np

# Read data
data_points = []
with open('step_02_sp_tau.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split()
        tau = float(parts[0])
        L = float(parts[3])
        data_points.append((tau, L))

# Filter out the point at τ = 0
data_points = [(tau, L) for tau, L in data_points if abs(tau) > 1e-10]

# Build design matrix using mpmath
n = len(data_points)
m = 4

A_mp = mp.matrix(n, m)
b_mp = mp.matrix(n, 1)

for i, (tau, L) in enumerate(data_points):
    tau_mp = mp.mpf(str(tau))
    L_mp = mp.mpf(str(L))
    
    A_mp[i, 0] = tau_mp**2
    A_mp[i, 1] = tau_mp**3
    A_mp[i, 2] = tau_mp**4
    A_mp[i, 3] = tau_mp**5
    b_mp[i, 0] = L_mp

ATA = A_mp.T * A_mp
ATb = A_mp.T * b_mp
coeffs_mp = mp.lu_solve(ATA, ATb)

# Compute expected values
expected_c2 = -mp.pi**2 / 9
expected_c3 = -16 * mp.zeta(3) / 9
expected_c4 = -37 * mp.pi**4 / 810
expected_c5 = -160 * mp.zeta(5) / 27

print("=" * 80)
print("COMPARISON: FITTED vs CLAIMED")
print("=" * 80)

def print_comparison(name, fitted, expected):
    diff = fitted - expected
    rel_err = abs(diff / expected)
    print(f"\n{name}:")
    print(f"  Fitted:   {mp.nstr(fitted, 25)}")
    print(f"  Expected: {mp.nstr(expected, 25)}")
    print(f"  Diff:     {mp.nstr(diff, 25)}")
    print(f"  Rel Err:  {mp.nstr(rel_err, 10)}")

print_comparison("c2 (τ²)", coeffs_mp[0], expected_c2)
print_comparison("c3 (τ³)", coeffs_mp[1], expected_c3)
print_comparison("c4 (τ⁴)", coeffs_mp[2], expected_c4)
print_comparison("c5 (τ⁵)", coeffs_mp[3], expected_c5)

# Compute residuals
print("\n" + "=" * 80)
print("RESIDUALS")
print("=" * 80)

residuals = []
for i, (tau, L) in enumerate(data_points):
    tau_mp = mp.mpf(str(tau))
    L_mp = mp.mpf(str(L))
    L_pred = sum(coeffs_mp[j] * tau_mp**(j+2) for j in range(4))
    residual = L_mp - L_pred
    residuals.append(float(residual))
    print(f"τ = {tau:+.3f}: {mp.nstr(residual, 8)}")

print(f"\nRMS residual: {np.sqrt(np.mean(np.array(residuals)**2)):.2e}")
print(f"Max |residual|: {np.max(np.abs(residuals)):.2e}")