import mpmath as mp
mp.mp.dps = 80
from mpmath import nsum, inf

def surface_area_monte_carlo(p, num_samples=2000000):
    """
    Estimate surface area using Monte Carlo integration.
    For testing/verification only.
    """
    import random
    mp.mp.dps = 80
    
    # Parameterize surface in first octant
    # Use uniform sampling on the surface
    
    total = mp.mpf(0)
    count = 0
    
    for _ in range(num_samples):
        # Random point in parameter space
        u = random.random()
        v = random.random()
        
        x = (u * v)**(1/p)
        y = (u * (1 - v))**(1/p)
        
        if x**p + y**p > 1:
            continue
        
        # Jacobian weight
        jac = (1/p**2) * u**(2/p - 1) * v**(1/p - 1) * (1 - v)**(1/p - 1)
        
        # Surface element contribution (approximate)
        # This is just for rough verification
        total += jac
        count += 1
    
    # Surface area estimate (very rough)
    area_estimate = 8 * total / count
    return area_estimate

def surface_area_direct(p, dps=50):
    """
    Direct integration using 2D quadrature with singularity handling.
    """
    old_dps = mp.mp.dps
    mp.mp.dps = dps + 20
    
    def f(v, u):
        """Integrand with u and v, u is the inner variable"""
        if u <= 0 or u >= 1 or v <= 0 or v >= 1:
            return mp.mpf(0)
        
        x = (u * v)**(1/mp.mpf(p))
        y = (u * (1 - v))**(1/mp.mpf(p))
        
        # Surface integrand: sqrt(1 + fx^2 + fy^2)
        fx2 = x**(2*(p-1)) * (1-u)**(2*(1-p)/p)
        fy2 = y**(2*(p-1)) * (1-u)**(2*(1-p)/p)
        
        surf = mp.sqrt(1 + fx2 + fy2)
        
        # Jacobian
        jac = (1/p**2) * u**(2/p - 1) * v**(1/p - 1) * (1 - v)**(1/p - 1)
        
        return surf * jac
    
    # Integrate using quadts (tanh-sinh) with error control
    # mp.quadts uses tanh-sinh which handles singularities well
    
    total = mp.quadts(lambda u: mp.quadts(lambda v: f(v, u), [0, 1]), [0, 1])
    
    mp.mp.dps = old_dps
    return 8 * total

# Quick test
print("Testing direct integration:")
for p in [1.5, 2.0]:
    mp.mp.dps = 60
    S = surface_area_direct(p, dps=50)
    mp.mp.dps = 60
    print(f"S({p}) = {S}")
    if p == 2.0:
        print(f"4*pi = {4*mp.pi}")