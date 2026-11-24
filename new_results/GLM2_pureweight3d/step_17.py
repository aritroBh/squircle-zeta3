import mpmath as mp
mp.mp.dps = 120

def surface_area_direct(p, dps=100):
    """Direct integration with higher precision."""
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

# Compute at many points symmetric around τ = 0
tau_values = []
for t in [0.005, 0.01, 0.015, 0.02, 0.025, 0.03]:
    tau_values.extend([t, -t])

print("Computing high-precision S(p) for fitting:")
print("=" * 80)

data_points = []
for tau in tau_values:
    p = 1/(1 - tau) if abs(tau) < 1 else mp.mpf(2)
    
    mp.mp.dps = 120
    S = surface_area_direct(p, dps=100)
    mp.mp.dps = 120
    
    S0 = 4 * mp.sqrt(3)
    factor = (2 - p) * S / S0
    L = mp.log(factor)
    
    data_points.append((float(tau), float(p), float(S), float(L)))
    print(f"τ = {tau:+.3f}: p = {float(p):.6f}, L = {float(L):.15e}")

print(f"\nComputed {len(data_points)} data points")

# Save
with open('step_03_high_precision.txt', 'w') as f:
    f.write("# tau, p, S(p), L(p)\n")
    for tau, p, S, L in data_points:
        f.write(f"{tau:.15e} {p:.15e} {S:.20e} {L:.20e}\n")

print("Saved to step_03_high_precision.txt")