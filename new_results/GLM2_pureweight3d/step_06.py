import mpmath as mp
mp.mp.dps = 80

def surface_area_optimized(p, dps=60):
    """
    Compute surface area S(p) of |x|^p + |y|^p + |z|^p = 1.
    Using coordinate transformation with proper singularity handling.
    """
    old_dps = mp.mp.dps
    mp.mp.dps = dps
    
    def integrand_uv(u, v):
        if u <= 0 or u >= 1 or v <= 0 or v >= 1:
            return mp.mpf(0)
        
        x = (u * v)**(1/p)
        y = (u * (1 - v))**(1/p)
        
        z_inner = 1 - u
        z_pow = z_inner**((1-p)/p)
        
        fx = -x**(p-1) * z_pow
        fy = -y**(p-1) * z_pow
        
        integrand = mp.sqrt(1 + fx**2 + fy**2)
        
        jacobian = (1/p**2) * u**(2/p - 1) * v**(1/p - 1) * (1 - v)**(1/p - 1)
        
        return integrand * jacobian
    
    def integral_v(u):
        return mp.quad(lambda v: integrand_uv(u, v), [0, 1])
    
    result = mp.quad(lambda u: integral_v(u), [0, 1])
    S_total = 8 * result
    
    mp.mp.dps = old_dps
    return S_total

# Save S(p) values to file
print("Computing S(p) for various p values:")
print("=" * 60)

results = []
d_values = [0.005, -0.005, 0.01, -0.01, 0.02, -0.02, 0.03, -0.03, 0.0]

for d in d_values:
    p = 1 + d
    mp.mp.dps = 80
    S = surface_area_optimized(p, dps=70)
    mp.mp.dps = 80
    
    if d == 0:
        expected = 4 * mp.sqrt(3)
        error = abs(S - expected) / expected
        print(f"\np = {p:.3f}:")
        print(f"  S(p) = {S}")
        print(f"  Expected = {expected}")
        print(f"  Relative error = {error}")
    else:
        print(f"\np = {p:.3f}:")
        print(f"  S(p) = {S}")
    
    results.append((p, S))

# Save to file
with open('sp_values.txt', 'w') as f:
    f.write("# p, S(p) computed to ~50 digits\n")
    for p, S in results:
        f.write(f"{p:.5f} {str(S)[:70]}\n")

print("\n" + "=" * 60)
print("Testing p = 2 (unit sphere):")
mp.mp.dps = 80
S2 = surface_area_optimized(2.0, dps=70)
mp.mp.dps = 80
print(f"S(2) = {S2}")
print(f"Expected 4*pi = {4*mp.pi}")
print(f"Relative error = {abs(S2 - 4*mp.pi)/(4*mp.pi)}")