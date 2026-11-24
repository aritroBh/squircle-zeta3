import mpmath as mp
mp.mp.dps = 120
import numpy as np

# Read high precision data
data_points = []
with open('step_03_high_precision.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split()
        tau = float(parts[0])
        L = float(parts[3])
        data_points.append((tau, L))

print(f"Loaded {len(data_points)} data points")

# Use mpmath for high precision fitting
mp.mp.dps = 120

n = len(data_points)
m = 6  # Fit τ, τ², τ³, τ⁴, τ⁵, τ⁶

A_mp = mp.matrix(n, m)
b_mp = mp.matrix(n, 1)

for i, (tau, L) in enumerate(data_points):
    tau_mp = mp.mpf(str(tau))
    L_mp = mp.mpf(str(L))
    
    # Full polynomial basis
    A_mp[i, 0] = tau_mp      # τ term (should be ~0)
    A_mp[i, 1] = tau_mp**2   # τ² term
    A_mp[i, 2] = tau_mp**3   # τ³ term
    A_mp[i, 3] = tau_mp**4   # τ⁴ term
    A_mp[i, 4] = tau_mp**5   # τ⁵ term
    A_mp[i, 5] = tau_mp**6   # τ⁶ term
    
    b_mp[i, 0] = L_mp

# Solve
ATA = A_mp.T * A_mp
ATb = A_mp.T * b_mp
coeffs_full = mp.lu_solve(ATA, ATb)

print("\n" + "=" * 80)
print("FULL POLYNOMIAL FIT")
print("=" * 80)
for j, name in enumerate(['τ', 'τ²', 'τ³', 'τ⁴', 'τ⁵', 'τ⁶']):
    print(f"c_{name} = {mp.nstr(coeffs_full[j], 30)}")

# Expected coefficients
expected = {
    'τ²': -mp.pi**2 / 9,
    'τ³': -16 * mp.zeta(3) / 9,
    'τ⁴': -37 * mp.pi**4 / 810,
    'τ⁵': -160 * mp.zeta(5) / 27,
}

print("\n" + "=" * 80)
print("COMPARISON WITH CLAIMED")
print("=" * 80)

for j, name in enumerate(['τ²', 'τ³', 'τ⁴', 'τ⁵']):
    idx = j + 1  # index in coeffs_full
    if name in expected:
        fitted = coeffs_full[idx]
        exp = expected[name]
        diff = fitted - exp
        rel_err = abs(diff / exp)
        print(f"\n{name}:")
        print(f"  Fitted:  {mp.nstr(fitted, 30)}")
        print(f"  Claimed: {mp.nstr(exp, 30)}")
        print(f"  Diff:    {mp.nstr(diff, 20)}")
        print(f"  Rel Err: {mp.nstr(rel_err, 10)}")

# Compute residuals with full polynomial
print("\n" + "=" * 80)
print("RESIDUALS WITH FULL POLYNOMIAL")
print("=" * 80)

residuals = []
for i, (tau, L) in enumerate(data_points):
    tau_mp = mp.mpf(str(tau))
    L_mp = mp.mpf(str(L))
    L_pred = sum(coeffs_full[j] * tau_mp**(j+1) for j in range(6))
    residual = L_mp - L_pred
    residuals.append(float(residual))

print(f"RMS residual: {np.sqrt(np.mean(np.array(residuals)**2)):.6e}")
print(f"Max |residual|: {np.max(np.abs(residuals)):.6e}")