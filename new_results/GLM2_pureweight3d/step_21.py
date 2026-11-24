import mpmath as mp
mp.mp.dps = 150

def surface_area_direct(p, dps=120):
    """Direct integration."""
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

# Compute many more points for robust fitting
# Use 20 points: τ = ±0.01, ±0.02, ..., ±0.10

tau_values = []
for t in [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05]:
    tau_values.extend([t, -t])

print(f"Computing {len(tau_values)} data points for robust fitting:")
print("=" * 80)

data_points = []
for tau in sorted(tau_values):
    p = 1/(1-mp.mpf(tau)) if abs(tau) < 0.999 else mp.mpf(2)
    
    mp.mp.dps = 150
    S = surface_area_direct(p, dps=120)
    mp.mp.dps = 150
    
    S0 = 4 * mp.sqrt(3)
    L = mp.log((2-p) * S / S0)
    
    data_points.append((float(tau), float(p), float(S), float(L)))
    print(f"τ = {tau:+.3f}: S = {float(S):.12f}, L = {float(L):.12e}")

# Save
with open('step_05_robust.txt', 'w') as f:
    f.write("# tau, p, S, L\n")
    for tau, p, S, L in data_points:
        f.write(f"{tau:.10e} {p:.15e} {S:.20e} {L:.20e}\n")

print(f"\nSaved {len(data_points)} points to step_05_robust.txt")