import mpmath as mp
mp.mp.dps = 150

def surface_area_direct(p, dps=120):
    """Direct integration with very high precision."""
    old_dps = mp.mp.dps
    mp.mp.dps = dps + 30
    
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

# Compute many points at very small τ values for accurate coefficient extraction
# Use Richardson extrapolation: compute at τ and τ/2, then extrapolate

print("Computing S(p) for Richardson extrapolation:")
print("=" * 80)

# Use τ = 0.02 and get coefficients from behavior there
tau_values = [0.02, 0.01, 0.005, 0.0025]

data = []
for tau in tau_values:
    p = 1/(1-mp.mpf(tau))
    
    mp.mp.dps = 150
    S = surface_area_direct(p, dps=120)
    mp.mp.dps = 150
    
    S0 = 4 * mp.sqrt(3)
    factor = (2 - p) * S / S0
    L = mp.log(factor)
    
    data.append((tau, float(p), float(S), float(L)))
    print(f"τ = {tau:.6f}: L = {float(L):.12e}")

# Save
with open('step_04_richardson.txt', 'w') as f:
    f.write("# tau, p, S, L (for Richardson extrapolation)\n")
    for tau, p, S, L in data:
        f.write(f"{tau:.10e} {p:.15e} {S:.20e} {L:.20e}\n")

print("\nSaved to step_04_richardson.txt")