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

# Filter out the point at τ = 0 (which gives L = 0 exactly)
data_points = [(tau, L) for tau, L in data_points if abs(tau) > 1e-10]

tau_arr = np.array([tau for tau, L in data_points])
L_arr = np.array([L for tau, L in data_points])

print(f"Using {len(data_points)} data points (τ ≠ 0)")
print()

# Build design matrix for polynomial: c2*τ² + c3*τ³ + c4*τ⁴ + c5*τ⁵
# Using high precision mpmath
mp.mp.dps = 100

# Convert to mpmath matrices
n = len(data_points)
m = 4  # 4 coefficients (τ², τ³, τ⁴, τ⁵)

# Build the design matrix A and vector b
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

# Solve using least squares: coeffs = (A^T A)^{-1} A^T b
ATA = A_mp.T * A_mp
ATb = A_mp.T * b_mp

# Solve the linear system
coeffs_mp = mp.lu_solve(ATA, ATb)

print("=" * 80)
print("Fitted coefficients (using mpmath least squares):")
print("=" * 80)
print(f"c2 (τ² coefficient) = {coeffs_mp[0]}")
print(f"c3 (τ³ coefficient) = {coeffs_mp[1]}")
print(f"c4 (τ⁴ coefficient) = {coeffs_mp[2]}")
print(f"c5 (τ⁵ coefficient) = {coeffs_mp[3]}")

# Compute expected values
expected_c2 = -mp.pi**2 / 9
expected_c3 = -16 * mp.zeta(3) / 9
expected_c4 = -37 * mp.pi**4 / 810
expected_c5 = -160 * mp.zeta(5) / 27

print("\n" + "=" * 80)
print("Expected coefficients (claimed formula):")
print("=" * 80)
print(f"c2 expected = {expected_c2}")
print(f"c3 expected = {expected_c3}")
print(f"c4 expected = {expected_c4}")
print(f"c5 expected = {expected_c5}")

# Compute relative errors
print("\n" + "=" * 80)
print("Relative errors (fitted vs expected):")
print("=" * 80)
err_c2 = abs((coeffs_mp[0] - expected_c2) / expected_c2)
err_c3 = abs((coeffs_mp[1] - expected_c3) / expected_c3)
err_c4 = abs((coeffs_mp[2] - expected_c4) / expected_c4)
err_c5 = abs((coeffs_mp[3] - expected_c5) / expected_c5)

print(f"Relative error in c2: {err_c2:.6e}")
print(f"Relative error in c3: {err_c3:.6e}")
print(f"Relative error in c4: {err_c4:.6e}")
print(f"Relative error in c5: {err_c5:.6e}")

# Compute residuals
residuals = []
for i, (tau, L) in enumerate(data_points):
    tau_mp = mp.mpf(str(tau))
    L_mp = mp.mpf(str(L))
    L_pred = coeffs_mp[0] * tau_mp**2 + coeffs_mp[1] * tau_mp**3 + \
             coeffs_mp[2] * tau_mp**4 + coeffs_mp[3] * tau_mp**5
    residual = L_mp - L_pred
    residuals.append(float(residual))

print("\n" + "=" * 80)
print("Residuals (L_observed - L_predicted):")
print("=" * 80)
for i, (tau, L) in enumerate(data_points):
    print(f"τ = {tau:+.3f}: residual = {residuals[i]:.2e}")

print(f"\nRMS residual: {np.sqrt(np.mean(np.array(residuals)**2)):.2e}")
print(f"Max absolute residual: {np.max(np.abs(residuals)):.2e}")