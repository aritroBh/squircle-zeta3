import mpmath as mp
mp.mp.dps = 100

def surface_area_fast(p, dps=80):
    """
    Compute surface area S(p) of |x|^p + |y|^p + |z|^p = 1.
    
    Fast implementation using singularity subtraction and transformation.
    """
    old_dps = mp.mp.dps
    mp.mp.dps = dps + 30
    
    p_mp = mp.mpf(p)
    inv_p = 1/p_mp
    
    # For the transformation: x^p = u*v, y^p = u*(1-v)
    # u in [0,1], v in [0,1]
    # This naturally handles the constraint x^p + y^p <= 1
    
    # The surface integrand is:
    # sqrt(1 + fx^2 + fy^2) * |∂(x,y)/∂(u,v)|
    # where z = (1 - u)^(1/p), and x = (u*v)^(1/p), y = (u*(1-v))^(1/p)
    
    # Compute ∂z/∂x and ∂z/∂y:
    # z^p = 1 - x^p - y^p
    # Taking partial derivative: p*z^(p-1) * ∂z/∂x = -p*x^(p-1)
    # So ∂z/∂x = -x^(p-1) * z^(1-p) = -x^(p-1) * (1-u)^((1-p)/p)
    
    # Near the boundary u=1: the surface flattens to z=0
    # For p>1: (1-p)/p < 0, so z^(1-p) diverges
    # Need to handle this singularity
    
    def integrand_plus_base(u, v):
        """Compute integrand minus singular part"""
        if u <= 0 or u >= 1 or v <= 0 or v >= 1:
            return mp.mpf(0)
        
        x = (u * v)**inv_p
        y = (u * (1 - v))**inv_p
        
        one_minus_u = 1 - u
        exp_factor = (1 - p_mp) * inv_p
        
        # fx = -x^(p-1) * (1-u)^((1-p)/p)
        # fy = -y^(p-1) * (1-u)^((1-p)/p)
        
        fx2 = x**(2*(p_mp - 1)) * one_minus_u**(2 * exp_factor)
        fy2 = y**(2*(p_mp - 1)) * one_minus_u**(2 * exp_factor)
        
        integrand = mp.sqrt(1 + fx2 + fy2)
        
        # Jacobian: |∂(x,y)/∂(u,v)| = (1/p^2) * u^(2/p-1) * v^(1/p-1) * (1-v)^(1/p-1)
        jacobian = inv_p**2 * u**(2*inv_p - 1) * v**(inv_p - 1) * (1 - v)**(inv_p - 1)
        
        # For singularity subtraction, compute the leading behavior near u=1
        # fx^2 + fy^2 ~ (x^(2(p-1)) + y^(2(p-1))) * (1-u)^(2*(1-p)/p)
        # = ((uv)^((p-1)/p) + (u(1-v))^((p-1)/p))^2 * (1-u)^(2*(1-p)/p)
        
        # Leading term: |∂z/∂x| (along x-direction) gives the edge contribution
        # For the edge where x^p = u, y^p = 0, the surface is z = (1-u)^(1/p)
        # Slope: |∂z/∂u| = (1/p) * (1-u)^(1/p - 1)
        
        # Actually, the singular contribution comes from sqrt(1 + fx^2 + fy^2)
        # when the gradient of z diverges at the edge u=1
        
        # For p > 1: exp_factor < 0, so near u=1:
        # sqrt(1 + fx^2 + fy^2) ~ sqrt(fx^2 + fy^2) ~ |∇z| ~ u^((p-1)/p) * (1-u)^(exp_factor)
        
        coeff = ((u*v)**((p_mp-1)*inv_p) + (u*(1-v))**((p_mp-1)*inv_p))
        
        # The singular part: sqrt(fx^2 + fy^2) approximately
        sing = coeff * one_minus_u**exp_factor
        
        return (integrand - sing) * jacobian
    
    def get_singular_contribution(p):
        """Compute analytic singular contribution"""
        # The integral of sqrt(fx^2 + fy^2) part
        # This equals ∫∫ u^((p-1)/p) * v^((p-1)/p) * (1-u)^((1-p)/p) * jacobian
        
        # For the singular part, the dominant contribution is:
        # ∫_{u=0}^1 ∫_{v=0}^1 (1-u)^((1-p)/p) * (1/p^2) * u^(2/p-1) * v^(1/p-1) * (1-v)^(1/p-1)
        
        # The coefficients from fx^2 and fy^2 create a weighting
        
        inv_p = 1/mp.mpf(p)
        
        # After careful analysis, the singular integral is:
        # ∫_{v=0}^1 ∫_{u=0}^1 sqrt(2) * u^((p-1)/p) * (v^((p-1)/p) + (1-v)^((p-1)/p))^^2 ...
        
        # This gets complicated. Let me use a numerical approximation.
        
        # Actually, the correct approach is to note that the integral of sqrt(fx^2+fy^2)
        # over the domain gives a finite contribution that can be computed by evaluating
        # the integral where we drop the "+1" in sqrt(1 + fx^2 + fy^2):
        
        def integrand_sing(u, v):
            if u <= 0 or u >= 1 or v <= 0 or v >= 1:
                return mp.mpf(0)
            
            coeff = ((u*v)**((p-1)/p) + (u*(1-v))**((p-1)/p))
            exp_factor = (1 - mp.mpf(p))/mp.mpf(p)
            
            sing = coeff * (1-u)**exp_factor
            jacobian = (1/p**2) * u**(2/p - 1) * v**(1/p - 1) * (1-v)**(1/p - 1)
            
            return sing * jacobian
        
        result = mp.quad(lambda u: mp.quad(lambda v: integrand_sing(u, v), [0, 1]), [0, 1])
        return result
    
    # Compute both parts
    integral_plus = mp.quad(lambda u: mp.quad(lambda v: integrand_plus_base(u, v), [0, 1]), [0, 1])
    sing_contribution = get_singular_contribution(p)
    
    # Need to add 1/sqrt contributions for points away from edge
    # Actually, I need to add back the integral of 1 (in the surface element)
    # which comes from the +1 in sqrt(1 + fx^2 + fy^2)
    
    # The "+1" contribution in the integrand:
    def integrand_base(u, v):
        if u <= 0 or u >= 1 or v <= 0 or v >= 1:
            return mp.mpf(0)
        jacobian = inv_p**2 * u**(2*inv_p - 1) * v**(inv_p - 1) * (1 - v)**(inv_p - 1)
        return jacobian
    
    integral_base = mp.quad(lambda u: mp.quad(lambda v: integrand_base(u, v), [0, 1]), [0, 1])
    
    # Total integral in first octant
    total_integral = integral_base + integral_plus + sing_contribution
    
    S_total = 8 * total_integral
    
    mp.mp.dps = old_dps
    return S_total

# Test
print("Testing fast implementation:")
for p in [2.0, 1.5, 1.2, 1.1]:
    mp.mp.dps = 80
    S = surface_area_fast(p, dps=70)
    mp.mp.dps = 80
    print(f"S({p}) = {S}")

print(f"\nExpected S(2) = 4*pi = {4*mp.pi}")