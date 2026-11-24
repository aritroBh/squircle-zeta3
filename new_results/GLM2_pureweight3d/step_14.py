import mpmath as mp
mp.mp.dps = 100

def surface_area_direct(p, dps=80):
    """Direct integration using 2D quadrature."""
    old_dps = mp.mp.dps
    mp.mp.dps = dps + 20
    
    p_mpf = mp.mpf(p)
    inv_p = 1/p_mpf
    
    def f(v, u):
        if u <= 0 or u >= 1 or v <= 0 or v >= 1:
            return mp.mpf(0)
        
        x = (u * v)**inv_p
        y = (u * (1 - v))**inv_p
        
        one_minus_u = 1 - u
        exp_factor = 2*(1-p_mpf)*inv_p
        
        fx2 = x**(2*(p_mpf-1)) * one_minus_u**exp_factor
        fy2 = y**(2*(p_mpf-1)) * one_minus_u**exp_factor
        
        surf = mp.sqrt(1 + fx2 + fy2)
        
        jac = inv_p**2 * u**(2*inv_p - 1) * v**(inv_p - 1) * (1 - v)**(inv_p - 1)
        
        return surf * jac
    
    total = mp.quadts(lambda u: mp.quadts(lambda v: f(v, u), [0, 1]), [0, 1])
    
    mp.mp.dps = old_dps
    return 8 * total

# Compute at finer grid near p = 1 to get more data points
# Use symmetric points in τ
tau_values = [0.0, 0.005, -0.005, 0.01, -0.01, 0.015, -0.015, 0.02, -0.02, 0.025, -0.025, 0.03, -0.03]

print("Computing S(p) for τ values:")
print("=" * 80)

data_points = []
for tau in tau_values:
    # τ = 1 - 1/p => p = 1/(1-τ)
    p = 1/(1 - tau) if abs(tau) < 1 else 2.0  # handle τ = 0
    
    mp.mp.dps = 100
    S = surface_area_direct(p, dps=85)
    mp.mp.dps = 100
    
    # Compute L(p) = log((2-p)*S/(4*sqrt(3)))
    if abs(p - 2) < 1e-10:
        # Special case for theoretical value
        L = mp.mpf(0)
    else:
        S0 = 4 * mp.sqrt(3)
        factor = (2 - p) * S / S0
        L = mp.log(factor)
    
    data_points.append((float(tau), float(p), float(S), float(L)))
    print(f"τ = {tau:+.3f}, p = {float(p):.6f}, S = {float(S):.12f}, L = {float(L):.15f}")

print(f"\nComputed {len(data_points)} data points")
print("Saving to file...")

with open('step_02_sp_tau.txt', 'w') as f:
    f.write("# tau, p, S(p), L(p) = log((2-p)*S/(4*sqrt(3)))\n")
    for tau, p, S, L in data_points:
        f.write(f"{tau:.10f} {p:.10f} {S:.18f} {L:.18f}\n")

print("Saved to step_02_sp_tau.txt")