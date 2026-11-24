import mpmath as mp
mp.mp.dps = 80

def surface_area_correct(p):
    """
    Compute surface area S(p) of |x|^p + |y|^p + |z|^p = 1.
    
    Surface integral: S = 8 * ∫∫_D sqrt(1 + (∂z/∂x)^2 + (∂z/∂y)^2) dx dy
    where D: x >= 0, y >= 0, x^p + y^p <= 1
    and z = (1 - x^p - y^p)^(1/p)
    
    Using x = (u*(1-t))^(1/p), y = (u*t)^(1/p)
    where u in [0, 1] and t in [0, 1/2] (by symmetry).
    """
    
    # First, compute the singular contribution analytically
    # As u -> 1 (edge of domain), the surface becomes flat
    # Contribution from singular part
    
    # Analyze: near u = 1, z ~ (1-u)^(1/p) -> 0
    # The singularity in the integrand comes from sqrt(1 + fx^2 + fy^2)
    # flattening to 1, but the Jacobian ~ (1-u)^(2/p - 2)
    
    # Actually, let me use a different parameterization
    # Use polar-like coordinates on the parametric surface
    
    def surface_parametric(theta, phi, p):
        """
        Parameterize the surface using:
        x = cos(theta)^(2/p) * sin(phi)^(2/p)
        y = sin(theta)^(2/p) * sin(phi)^(2/p)  
        z = cos(phi)^(2/p)
        normalized appropriately.
        """
        # Not quite right either...
        pass
    
    # Let's use the direct formula for surface area in first octant
    # S_1 = ∫∫_D sqrt(1 + fx^2 + fy^2) dx dy
    
    def integrand_double(x, y, p):
        """Compute integrand sqrt(1 + fx^2 + fy^2) at point (x, y)"""
        if x <= 0 or y <= 0:
            return mp.mpf(1)  # boundary
        
        xp = x**p
        yp = y**p
        
        if xp + yp >= 1:
            return mp.mpf(0)  # outside domain
        
        # z = (1 - x^p - y^p)^(1/p)
        # ∂z/∂x = -x^(p-1) * z^(1-p) = -x^(p-1) * (1 - x^p - y^p)^((1-p)/p)
        # ∂z/∂y = -y^(p-1) * z^(1-p)
        
        z_inner = 1 - xp - yp
        z_pow = z_inner**((1-p)/p)
        
        fx = -x**(p-1) * z_pow
        fy = -y**(p-1) * z_pow
        
        result = mp.sqrt(1 + fx**2 + fy**2)
        return result
    
    # Now integrate. Use transformation to handle singularity at x^p + y^p = 1
    # Let u = x^p + y^p, v = x^p / (x^p + y^p) (when u > 0)
    # Then x = (u*v)^(1/p), y = (u*(1-v))^(1/p)
    # dx dy = (1/p^2) * u^(2/p - 1) * v^(1/p - 1) * (1-v)^(1/p - 1) du dv
    
    def integrand_uv(u, v, p):
        """Transformed integrand in (u, v) coordinates"""
        if u <= 0 or u >= 1 or v <= 0 or v >= 1:
            return mp.mpf(0)
        
        x = (u * v)**(1/p)
        y = (u * (1 - v))**(1/p)
        
        # Original integrand
        inner = integrand_double(x, y, p)
        
        # Jacobian
        jacobian = (1/p**2) * u**(2/p - 1) * v**(1/p - 1) * (1 - v)**(1/p - 1)
        
        return inner * jacobian
    
    # Handle singularity at u = 1
    # Near u = 1, z ~ (1-u)^(1/p) -> 0
    # fx, fy -> finite values
    # integrand ~ sqrt(1 + const*(1-u)^(2*(1-p)/p))
    # For p > 1, this diverges but is integrable
    
    # Subtract and add the leading singular behavior
    # sqrt(1 + c*(1-u)^alpha) ~ 1 + (c/2)*(1-u)^alpha for small (1-u)
    
    def integrand_uv_subtracted(u, v, p):
        """Integrand with singularity subtraction"""
        if u <= 0 or u >= 1 or v <= 0 or v >= 1:
            return mp.mpf(0)
        
        x = (u * v)**(1/p)
        y = (u * (1 - v))**(1/p)
        
        z_inner = 1 - u
        z_pow = z_inner**((1-p)/p)
        
        fx = -x**(p-1) * z_pow
        fy = -y**(p-1) * z_pow
        
        inner = mp.sqrt(1 + fx**2 + fy**2)
        
        # Singular contribution (to be added back analytically)
        # As u -> 1: fx^2 + fy^2 ~ x^(2p-2) + y^(2p-2) times (1-u)^(2(1-p)/p)
        # For the singularity subtraction, use: 1/sqrt(1-u) * (const based on p)
        
        sing_factor = (x**(p-1))**2 + (y**(p-1))**2
        exponent = 2*(1-p)/p
        
        # For p > 1: exponent < 0, so singularity as u -> 1
        if p > 1:
            # Near u=1: fx^2 + fy^2 ~ sing_factor * (1-u)^exp where exp < 0
            # sqrt(1 + fx^2 + fy^2) ~ sqrt(fx^2 + fy^2) ~ sqrt(sing_factor) * (1-u)^(exp/2)
            # This is the leading singular behavior
            sing = mp.sqrt(sing_factor) * (1 - u)**(exponent/2)
        else:
            sing = mp.mpf(1)
        
        jacobian = (1/p**2) * u**(2/p - 1) * v**(1/p - 1) * (1 - v)**(1/p - 1)
        
        return (inner - sing) * jacobian
    
    # Actually, let me just integrate directly with careful handling
    
    # Direct integration over triangle: x in [0, 1], y in [0, (1-x^p)^(1/p)]
    def inner_y(x, p):
        y_max = (1 - x**p)**(1/p)
        def f(y):
            return integrand_double(x, y, p)
        return mp.quad(f, [0, y_max])
    
    def inner_x(p):
        return mp.quad(lambda x: inner_y(x, p), [0, 1])
    
    S_first_octant = inner_x(p)
    S_total = 8 * S_first_octant
    
    return S_total

# Test
print("Testing surface area calculation:")
S2 = surface_area_correct(2)
print(f"S(2) = {S2}")
print(f"Expected 4*pi = {4*mp.pi}")