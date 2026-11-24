import mpmath as mp
mp.mp.dps = 80

def surface_area_optimized(p, dps=50):
    """
    Compute surface area S(p) of |x|^p + |y|^p + |z|^p = 1.
    
    Using coordinate transformation to handle edge singularity.
    """
    old_dps = mp.mp.dps
    mp.mp.dps = dps
    
    # Use transformation: x = (u*v)^(1/p), y = (u*(1-v))^(1/p)
    # u in [0, 1], v in [0, 1]
    # Jacobian: J = (1/p^2) * u^(2/p - 1) * v^(1/p - 1) * (1-v)^(1/p - 1)
    
    def integrand_uv(u, v):
        """Transformed integrand in (u, v) coordinates"""
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
    
    # Split integration to handle singularities:
    # u near 0: u^(2/p - 1) singular for p > 2
    # u near 1: edge of surface, z = 0
    # v near 0 or 1: x^p or y^p contribution
    
    # Use tanh-sinh which handles endpoint singularities well
    
    def integral_v(u):
        return mp.quad(lambda v: integrand_uv(u, v), [0, 1])
    
    result = mp.quad(lambda u: integral_v(u), [0, 1])
    S_total = 8 * result
    
    mp.mp.dps = old_dps
    return S_total

# Test with p = 2 and p close to 1
print("Testing optimized surface area calculation:")
for p_test in [2.0, 1.5]:
    mp.mp.dps = 80
    S_test = surface_area_optimized(p_test, dps=70)
    mp.mp.dps = 80
    print(f"S({p_test}) = {S_test}")

print(f"\nExpected S(2) = 4*pi = {4*mp.pi}")
print(f"Expected S(1) = 4*sqrt(3) = {4*mp.sqrt(3)}")