import mpmath as mp
mp.mp.dps = 80
import numpy as np
from numpy.polynomial import polynomial as P

# Use the results already computed
results_data = """
1.005 6.962827472651343110847820492729664639603535281138012392327062728794635322047557058990716488679255571
0.995 6.893545510163452955986969972011966965763878758906176051180319816870103438244348728721135971074305679
1.010 6.997417989197152598961926014827674218735387989531841880859762847407137000900279984892250478860521261
0.990 6.858854480630419957393748911418098463824159640415041889698569347058257877662163988664326301511454716
1.020 7.066496490268251273907312209703136334581962230143539459609640371599784251188675294097641749599212251
0.980 6.789372730228341927674027348153347515965045205727414911142650897372367076306631695725270947604616769
1.030 7.135435383008283592517837599461182350016922833680139833466543149839800813722176329335572374525988870
0.970 6.719757492661952332608784377782465815057225885943284049550171686842912576254707577116833694993540153
1.000 6.928203230275509174109785366023489467771221015241522512223227917807732067635200148324584747028994303
"""

# Parse results
results = {}
for line in results_data.strip().split('\n'):
    parts = line.split()
    p = mp.mpf(parts[0])
    S = mp.mpf(parts[1])
    results[float(p)] = float(S)

# Create arrays for fitting
p_arr = np.array(sorted(results.keys()))
S_arr = np.array([results[p] for p in p_arr])

# Compute τ = 1 - 1/p
tau_arr = 1 - 1/p_arr

# Compute L(p) = log((2-p)*S(p)/(4*sqrt(3)))
S0 = 4 * mp.sqrt(3)  # S(1)
L_arr = np.log((2 - p_arr) * S_arr / float(S0))

print("Data for fitting:")
print("p          τ          L(p)")
print("-" * 50)
for i in range(len(p_arr)):
    print(f"{p_arr[i]:.5f}  {tau_arr[i]:.5f}  {L_arr[i]:.12f}")

# Fit polynomial L = c2*τ^2 + c3*τ^3 + c4*τ^4 + c5*τ^5
# Using numpy polyfit with higher precision
degree = 5

# Prepare design matrix
X = np.zeros((len(tau_arr), degree))
for i, tau in enumerate(tau_arr):
    for j in range(degree):
        X[i, j] = tau**(j+1)  # Note: j+1 gives τ, τ^2, ..., but we want τ^2, τ^3, ...
# Actually, let me match the polynomial: c2*τ^2 + c3*τ^3 + c4*τ^4 + c5*τ^5
# So coefficients are: 0 + 0*τ + c2*τ^2 + c3*τ^3 + c4*τ^4 + c5*τ^5

X_correct = np.zeros((len(tau_arr), degree))
for i, tau in enumerate(tau_arr):
    for j in range(degree):
        X_correct[i, j] = tau**(j+2)  # τ^2, τ^3, τ^4, τ^5, τ^6

# Fit using least squares
coeffs, residuals, rank, s = np.linalg.lstsq(X_correct, L_arr, rcond=None)

print("\n" + "=" * 50)
print("Fitted coefficients:")
print(f"c2 (τ² coefficient) = {coeffs[0]:.15f}")
print(f"c3 (τ³ coefficient) = {coeffs[1]:.15f}")
print(f"c4 (τ⁴ coefficient) = {coeffs[2]:.15f}")
print(f"c5 (τ⁵ coefficient) = {coeffs[3]:.15f}")

# Expected values
expected_c2 = -mp.pi**2 / 9
expected_c3 = -16 * mp.zeta(3) / 9
expected_c4 = -37 * mp.pi**4 / 810
expected_c5 = -160 * mp.zeta(5) / 27

print("\n" + "=" * 50)
print("Expected coefficients (claimed formula):")
print(f"c2 expected = {float(expected_c2):.15f}")
print(f"c3 expected = {float(expected_c3):.15f}")
print(f"c4 expected = {float(expected_c4):.15f}")
print(f"c5 expected = {float(expected_c5):.15f}")