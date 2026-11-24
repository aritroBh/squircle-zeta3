import mpmath as mp
mp.mp.dps = 80

def surface_area_correct(p):
    """
    Compute surface area S(p) of |x|^p + |y|^p + |z|^p = 1.
    
    Surface integral: S = 8 * ∫∫_D sqrt(1 + (∂z/∂x)^2 + (∂z/∂y)^2) dx dy
    where D: x >= 0, y >= 0, x^p + y^p <= 1
    and z = (1 - x^p - y^p)^(1/p)
    """
    
    def integrand_double(x, y, p):
        """Compute integrand sqrt(1 + fx^2 + fy^2) at point (x, y)"""
        if x <= 0 or y <= 0:
            return mp.mpf(1)
        
        xp = x**p
        yp = y**p
        
        if xp + yp >= 1:
            return mp.mpf(0)
        
        z_inner = 1 - xp - yp
        z_pow = z_inner**((1-p)/p)
        
        fx = -x**(p-1) * z_pow
        fy = -y**(p-1) * z_pow
        
        result = mp.sqrt(1 + fx**2 + fy**2)
        return result
    
    def inner_y(x, p):
        y_max = (1 - x**p)**(1/p)
        return mp.quad(lambda y: integrand_double(x, y, p), [0, y_max])
    
    def inner_x(p):
        return mp.quad(lambda x: inner_y(x, p), [0, 1])
    
    S_first_octant = inner_x(p)
    S_total = 8 * S_first_octant
    
    return S_total

# Tests
print("Testing surface area calculation:")
S2 = surface_area_correct(2)
print(f"S(2) = {S2}")
print(f"Expected 4*pi = {4*mp.pi}")
print(f"Relative error: {abs(S2 - 4*mp.pi)/(4*mp.pi)}")
print()

# Test p = 1 (octahedron): S(1) should be 4*sqrt(3)
for p_test in [1.01, 1.1]:
    S_test = surface_area_correct(p_test)
    print(f"S({p_test}) = {S_test}")

print("\nExpected S(1) = 4*sqrt(3) =", 4*mp.sqrt(3))