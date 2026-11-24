import mpmath as mp
mp.mp.dps = 100

def surface_area_hybrid(p, dps=80):
    """
    Compute surface area S(p) of |x|^p + |y|^p + |z|^p = 1.
    
    Uses integration over first octant:
    S = 8 * ∫∫_D sqrt(1 + fx^2 + fy^2) dx dy
    
    where z = (1 - x^p - y^p)^(1/p)
    """
    old_dps = mp.mp.dps
    mp.mp.dps = dps + 20
    
    one_p = 1/mp.mpf(p)
    two_p_minus_1 = 2*one_p - 1
    one_p_minus_1 = one_p - 1
    
    if p >= 1:
        # Use u = x^p + y^p transformation
        def integrand_uv(u, v):
            if u <= 0 or u >= 1 or v <= 0 or v >= 1:
                return mp.mpf(0)
            
            x = (u * v)**one_p
            y = (u * (1 - v))**one_p
            
            z_inner = 1 - u
            z_pow = z_inner**((1-p)*one_p)
            
            fx = -x**(p-1) * z_pow
            fy = -y**(p-1) * z_pow
            
            integrand = mp.sqrt(1 + fx**2 + fy**2)
            
            jacobian = (one_p**2) * u**two_p_minus_1 * v**one_p_minus_1 * (1 - v)**one_p_minus_1
            
            return integrand * jacobian
    else:
        # For p < 1, use direct integration with exponent transformation
        # x = u, y = v*u where u in [0,1], v in [0, 1]
        # Jacobian: u
        
        def integrand_uv(u, v):
            if u <= 0 or u >= 1:
                return mp.mpf(0)
            
            x = u**(1/p) * v**(1/p)
            y = u**(1/p) * (1-v)**(1/p)
            
            # Actually, let me use: x = u^(1/p), y = v * (1 - x^p)^(1/p) = v * (1-u)^(1/p)
            # Jacobian = (1-u)^(1/p - 1)
            pass
        # This is getting complicated. Let me use a simpler but slower approach for p < 1
        
        def integrand_xy(x, y):
            if x <= 0 or y <= 0:
                return mp.mpf(0)
            if x**p + y**p >= 1:
                return mp.mpf(0)
            
            z_inner = 1 - x**p - y**p
            z_pow = z_inner**((1-p)*one_p)
            
            fx = -x**(p-1) * z_pow
            fy = -y**(p-1) * z_pow
            
            return mp.sqrt(1 + fx**2 + fy**2)
        
        def inner_integral(x):
            y_max = (1 - x**p)**one_p
            return mp.quad(lambda y: integrand_xy(x, y), [0, y_max])
        
        result = mp.quad(inner_integral, [0, 1])
        S_total = 8 * result
        mp.mp.dps = old_dps
        return S_total
    
    def integral_v(u):
        return mp.quad(lambda v: integrand_uv(u, v), [0, 1])
    
    result = mp.quad(lambda u: integral_v(u), [0, 1])
    S_total = 8 * result
    
    mp.mp.dps = old_dps
    return S_total

# Test with various p values
print("Testing hybrid implementation:")
test_p_values = [2.0, 1.5, 1.1, 1.05, 1.01, 1.0, 0.99, 0.95, 0.9]

results = []
for p in test_p_values:
    mp.mp.dps = 100
    S = surface_area_hybrid(p, dps=80)
    mp.mp.dps = 100
    results.append((p, S))
    
    if p == 2.0:
        expected = 4 * mp.pi
    elif p == 1.0:
        expected = 4 * mp.sqrt(3)
    else:
        expected = None
    
    if expected:
        error = abs(S - expected) / expected
        print(f"p = {p:.2f}: S(p) = {S}")
        print(f"         Expected = {expected}")
        print(f"         Relative error = {error:.2e}")
    else:
        print(f"p = {p:.2f}: S(p) = {S}")

# Save results
with open('sp_values_hybrid.txt', 'w') as f:
    for p, S in results:
        f.write(f"{p:.10f} {str(S)[:90]}\n")