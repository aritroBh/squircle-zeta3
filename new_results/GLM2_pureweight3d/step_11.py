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
        
        # Surface integrand: sqrt(1 + fx^2 + fy^2)
        one_minus_u = 1 - u
        exp_factor = 2*(1-p_mpf)*inv_p
        
        fx2 = x**(2*(p_mpf-1)) * one_minus_u**exp_factor
        fy2 = y**(2*(p_mpf-1)) * one_minus_u**exp_factor
        
        surf = mp.sqrt(1 + fx2 + fy2)
        
        # Jacobian
        jac = inv_p**2 * u**(2*inv_p - 1) * v**(inv_p - 1) * (1 - v)**(inv_p - 1)
        
        return surf * jac
    
    total = mp.quadts(lambda u: mp.quadts(lambda v: f(v, u), [0, 1]), [0, 1])
    
    mp.mp.dps = old_dps
    return 8 * total

# Compute S(p) for required d values
print("Computing S(p) for p = 1 + d:")
print("=" * 70)

d_values = [0.005, -0.005, 0.01, -0.01, 0.02, -0.02, 0.03, -0.03, 0.0]
results = {}

for d in d_values:
    p = 1 + d
    print(f"\nComputing S({p:.3f})...")
    mp.mp.dps = 100
    S = surface_area_direct(p, dps=80)
    mp.mp.dps = 100
    
    results[p] = S
    print(f"S({p:.3f}) = {S}")
    
    if p == 2.0:
        print(f"Expected 4*pi = {4*mp.pi}")
        print(f"Rel error = {abs(S - 4*mp.pi)/(4*mp.pi):.2e}")
    elif abs(p - 1.0) < 1e-10:
        print(f"Expected 4*sqrt(3) = {4*mp.sqrt(3)}")
        print(f"Rel error = {abs(S - 4*mp.sqrt(3))/(4*mp.sqrt(3)):.2e}")

# Save to file
print("\n" + "=" * 70)
print("Saving results to file...")

with open('step_01_sp_values.txt', 'w') as f:
    f.write("# Computed surface areas S(p) for |x|^p + |y|^p + |z|^p = 1\n")
    f.write("# Using direct 2D quadrature (tanh-sinh)\n")
    f.write("# Format: p, S(p)\n")
    for p in sorted(results.keys()):
        S = results[p]
        f.write(f"{p:.10f} {str(S)[:85]}\n")

print("Done!")