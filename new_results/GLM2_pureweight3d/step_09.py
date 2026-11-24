import mpmath as mp
mp.mp.dps = 80

def surface_area_fast2(p, dps=60):
    """
    Compute surface area S(p) of |x|^p + |y|^p + |z|^p = 1.
    
    Uses integration over first octant with proper transformations.
    """
    old_dps = mp.mp.dps
    mp.mp.dps = dps + 10
    
    p_mp = mp.mpf(p)
    inv_p = 1/p_mp
    
    def inner_integrand(u, v):
        """
        Integrand in (u, v) coordinates where:
        x^p = u*v, y^p = u*(1-v)
        u in [0,1], v in [0,1]
        """
        if u <= 0 or u >= 1 or v <= 0 or v >= 1:
            return mp.mpf(0)
        
        x = (u * v)**inv_p
        y = (u * (1 - v))**inv_p
        
        # Compute the surface element
        # z = (1 - u)^(1/p)
        # ∂z/∂x = -x^(p-1) * (1-u)^((1-p)/p)
        # ∂z/∂y = -y^(p-1) * (1-u)^((1-p)/p)
        
        one_minus_u = 1 - u
        exp_factor = (1 - p_mp) * inv_p
        
        grad_z_sq = (x**(p_mp-1) + y**(p_mp-1))**2 * one_minus_u**(2*exp_factor)
        
        # Wait, this should be x^(2(p-1)) + y^(2(p-1)), not the square of sum
        fx2 = x**(2*(p_mp-1)) * one_minus_u**(2*exp_factor)
        fy2 = y**(2*(p_mp-1)) * one_minus_u**(2*exp_factor)
        
        surface_element = mp.sqrt(1 + fx2 + fy2)
        
        # Jacobian of transformation
        # x = (uv)^(1/p), y = (u(1-v))^(1/p)
        # ∂x/∂u = (1/p) * v^(1/p) * u^(1/p - 1)
        # ∂x/∂v = (1/p) * u^(1/p) * v^(1/p - 1)
        # ∂y/∂u = (1/p) * (1-v)^(1/p) * u^(1/p - 1)
        # ∂y/∂v = -(1/p) * u^(1/p) * (1-v)^(1/p - 1)
        # |det| = (1/p^2) * u^(2/p - 1) * (v^(1/p-1) * (1-v)^(1/p) + v^(1/p) * (1-v)^(1/p-1))
        #       = (1/p^2) * u^(2/p - 1) * v^(1/p-1) * (1-v)^(1/p-1)
        
        jacobian = inv_p**2 * u**(2*inv_p - 1) * v**(inv_p - 1) * (1 - v)**(inv_p - 1)
        
        return surface_element * jacobian
    
    # Use mp.quadts for tanh-sinh integration
    def inner_v(u):
        return mp.quad(lambda v: inner_integrand(u, v), [0, 1])
    
    result = mp.quad(lambda u: inner_v(u), [0, 1])
    
    S_total = 8 * result
    
    mp.mp.dps = old_dps
    return S_total

# Test
print("Testing fast2 implementation:")
for p in [2.0, 1.5, 1.2, 1.05]:
    mp.mp.dps = 80
    S = surface_area_fast2(p, dps=70)
    mp.mp.dps = 80
    if p == 2.0:
        print(f"S({p}) = {S}")
        print(f"Expected 4*pi = {4*mp.pi}")
        print(f"Error = {abs(S - 4*mp.pi)/(4*mp.pi):.2e}")
    else:
        print(f"S({p}) = {S}")

print(f"\nExpected S(1) = 4*sqrt(3) = {4*mp.sqrt(3)}")