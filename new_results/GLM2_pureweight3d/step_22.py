import mpmath as mp
mp.mp.dps = 100

# Summary of computed values from earlier runs
print("=" * 80)
print("TASK GLM2: INDEPENDENT CHECK OF PURE WEIGHT CLAIM IN 3D")
print("=" * 80)

# Computed S(p) values from step_03_high_precision.txt
data_str = """
τ = +0.005: L = -2.768549181871190e-05
τ = -0.005: L = -2.715120589737359e-05
τ = +0.010: L = -1.118443783270199e-04
τ = -0.010: L = -1.075691687521279e-04
τ = +0.015: L = -2.541824593684432e-04
τ = -0.015: L = -2.397484352462582e-04
τ = +0.020: L = -4.564770729905356e-04
τ = -0.020: L = -4.222458388892841e-04
τ = +0.025: L = -7.205796312082933e-04
τ = -0.025: L = -6.536783590352181e-04
τ = +0.030: L = -1.048418268781002e-03
τ = -0.030: L = -9.327210322811750e-04
"""

# Fitted coefficients from earlier (full polynomial with RMS 6e-12)
print("\nFITTED COEFFICIENTS (from least squares on 12 data points):")
print("c_τ  ≈ -1.84e-9 (should be zero)")
print("c_τ² = -1.096622719(8)")
print("c_τ³ = -2.13697248(9)")
print("c_τ⁴ = -4.4494933(58)")
print("c_τ⁵ = -6.185340(41)")

# Expected values
expected_c2 = -mp.pi**2 / 9
expected_c3 = -16 * mp.zeta(3) / 9
expected_c4 = -37 * mp.pi**4 / 810
expected_c5 = -160 * mp.zeta(5) / 27

print("\nEXPECTED COEFFICIENTS (claimed):")
print(f"c_τ² = {mp.nstr(expected_c2, 30)}")
print(f"c_τ³ = {mp.nstr(expected_c3, 30)}")
print(f"c_τ⁴ = {mp.nstr(expected_c4, 30)}")
print(f"c_τ⁵ = {mp.nstr(expected_c5, 30)}")

# Numerical comparisons
print("\nRELATIVE ERRORS:")
print("τ²: ~ 7.3e-9")
print("τ³: ~ 8.2e-6") 
print("τ⁴: ~ 1.3e-5")
print("τ⁵: ~ 6.6e-3")

print("\nRMS residual: 6.25e-12")
print("Max |residual|: 8.87e-12")

# Verification cases
print("\n" + "=" * 80)
print("VERIFICATION CASES:")
print("=" * 80)
print(f"S(1) = 6.9282032302755091... (matches 4√3 to 70+ digits)")
print(f"S(2) = 12.5663706143591729... (matches 4π to 38+ digits)")

print("\n" + "=" * 80)
print("FILES CREATED:")
print("=" * 80)
print("step_01_sp_values.txt - S(p) for p = 1+d")
print("step_02_sp_tau.txt - L(τ) values")
print("step_03_high_precision.txt - High precision data")
print("step_04_richardson.txt - Richardson extrapolation data")
print("step_05_robust.txt - (partial) Additional data")