import mpmath as mp
mp.mp.dps = 100

def surface_area_direct(p, dps=80):
    """
    Direct integration using 2D quadrature with singularity handling.
    """
    old_dps = mp.mp.dps
    mp.mp.dps = dps + 20
    
    p_mpf = mp.mpf(p)
    inv_p = 1/p_mpf
    
    def f(v, u):
        """Integrand with u and v swapped for proper integration order"""
        if u <= 0 or u >= 1 or v <= 0 or v >= 1:
            return mp.mpf(0)
        
        x = (u * v)**inv_p
        y = (u * (1 - v))**inv_p
        
        one_minus_u = 1 - u
        exp_factor = 2*(1-p_mpf)*inv_p
        
        fx2 = x**(2*(p_mpf-1)) * one_minus_u**exp_factor
        fy2 = y**(2*(p_mpf-1)) * one_minus_u**exp_factor
        
        surf = mp.sqrt(1 + fx2 + fy2)
        
        jac = inv_p**2 * u**(2*inv_p - 1) * v**(inv_p - 1) * (1 - v)**(inv_p - 1)
        
        return surf * jac
    
    total = mp.quadts(lambda u: mp.quadts(lambda v: f(v, u), [0, 1]), [0, 1])
    
    mp.mp.dps = old_dps
    return 8 * total

# Compute S(p) for required d values
d_values = [0.005, -0.005, 0.01, -0.01, 0.02, -0.02, 0.03, -0.03]
results = {}

print("Computing S(p) for p = 1 + d:")
print("=" * 70)

for d in d_values:
    p = 1 + d
    mp.mp.dps = 100
    S = surface_area_direct(p, dps=80)
    mp.mp.dps = 100
    
    results[p] = S
    print(f"S({p:+.3f}) = {S}")

# Verify special cases
print("\n" + "=" * 70)
print("Verification:")
mp.mp.dps = 100
S1 = surface_area_direct(1.0, dps=80)
S2 = surface_area_direct(2.0, dps=80)
mp.mp.dps = 100

print(f"S(1) = {S1}")
print(f"Expected 4*sqrt(3) = {4*mp.sqrt(3)}")
err1 = abs(S1 - 4*mp.sqrt(3))/(4*mp.sqrt(3))
print(f"Rel error = {mp.nstr(err1, 10)}")

print(f"\nS(2) = {S2}")
print(f"Expected 4*pi = {4*mp.pi}")
err2 = abs(S2 - 4*mp.pi)/(4*mp.pi)
print(f"Rel error = {mp.nstr(err2, 10)}")

# Save to file
with open('step_01_sp_values.txt', 'w') as f:
    f.write("# Computed surface areas S(p)\n")
    for p in sorted(results.keys()):
        S = results[p]
        f.write(f"{p:.10f} {str(S)[:85]}\n")

print("\nResults saved to step_01_sp_values.txt")