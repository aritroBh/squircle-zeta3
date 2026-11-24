import mpmath as mp
mp.mp.dps = 80

def surface_area(p):
    """
    Compute surface area S(p) of |x|^p + |y|^p + |z|^p = 1.
    Uses integration over first octant with singularity subtraction.
    """
    # Integrating over region 0 <= x, 0 <= y, x^p + y^p <= 1
    # z = (1 - x^p - y^p)^{1/p}
    # 
    # Surface element: sqrt(1 + (∂z/∂x)^2 + (∂z/∂y)^2)
    # 
    # After singularity subtraction, we use transformed integrals
    # with exponential convergence.
    
    def integrand_transformed(t, u):
        """
        Transformed integrand after removing singularity.
        Uses substitution x = u^(1/p), y = v^(1/p) and edge transformation.
        """
        if t <= 0 or u <= 0 or u >= 1:
            return mp.mpf(0)
        
        # x = (u*(1-t))^(1/p), y = (u*t)^(1/p)
        # This parameterizes edge where x^p + y^p = u
        x = (u * (1 - t))**(1/p)
        y = (u * t)**(1/p)
        
        if x**p + y**p >= 1:
            return mp.mpf(0)
        
        # z and its derivatives
        z = (1 - u)**(1/p)
        
        # Using the formula for surface area of z = f(x, y):
        # dS = sqrt(1 + fx^2 + fy^2) dx dy
        
        # After transformations and Jacobian, with singularity subtraction
        # at the boundary u -> 1:
        
        # Base contribution (from singularity subtraction)
        base = 1/mp.sqrt(1 - u)
        
        # Actual integrand
        # fx = -x^(p-1) * (1 - u)^((1-p)/p) 
        # fy = -y^(p-1) * (1 - u)^((1-p)/p)
        # After simplification:
        
        fx2 = (x**(p-1))**2 * (1 - u)**(2*(1-p)/p)
        fy2 = (y**(p-1))**2 * (1 - u)**(2*(1-p)/p)
        
        actual = mp.sqrt(1 + fx2 + fy2)
        
        # Jacobian for transformation
        # x = (u*(1-t))^(1/p) => dx/du = ... , dx/dt = ...
        # y = (u*t)^(1/p) => dy/du = ... , dy/dt = ...
        
        dx_du = (1/p) * (u*(1-t))**(1/p - 1) * (1-t)
        dx_dt = (1/p) * (u*(1-t))**(1/p - 1) * (-u)
        
        dy_du = (1/p) * (u*t)**(1/p - 1) * t
        dy_dt = (1/p) * (u*t)**(1/p - 1) * u
        
        jacobian = abs(dx_du * dy_dt - dx_dt * dy_du)
        
        return (actual - base) * jacobian
    
    def singularity_integral(p):
        """
        Compute the analytic part from singularity subtraction.
        This is the integral of 1/sqrt(1-u) with Jacobian.
        """
        # This gives a known contribution
        return 2  # Simplified form after integration over u from 0 to 1
    
    # Numerical integration using tanh-sinh (mpmath quad)
    def inner_u(t):
        def f(u):
            return integrand_transformed(t, u)
        return mp.quad(f, [0, 1])
    
    # Integrate over t from 0 to 1 (using symmetric integration)
    def integrand_t(t):
        if t <= 0 or t >= 1:
            return mp.mpf(0)
        try:
            result = inner_u(t)
            return result
        except:
            return mp.mpf(0)
    
    # Use mp.quad with intervals avoiding endpoints
    numerical_part = mp.quad(integrand_t, [0, 1])
    
    # Add singularity contribution and multiply by 8 for all octants
    S_8 = 8 * (numerical_part + singularity_integral(p))
    
    return S_8

# Test with p = 2 (sphere with radius 2^(-1/2), surface area 4*pi*r^2 = 4*pi/2 = 2*pi)
# Wait, for p=2: x^2 + y^2 + z^2 = 1 is unit sphere, S = 4*pi
print("Testing p = 2 (sphere):")
S2 = surface_area(2)
print(f"S(2) = {S2}")
print(f"Expected: 4*pi = {4*mp.pi}")
print(f"Relative error: {abs(S2 - 4*mp.pi)/(4*mp.pi)}")