import mpmath as mp
mp.mp.dps = 150

# Read data
data = []
with open('step_04_richardson.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split()
        tau = float(parts[0])
        L = float(parts[3])
        data.append((tau, L))

print("Richardson extrapolation for L(τ):")
print("=" * 80)

# L(τ) = c2*τ² + c3*τ³ + c4*τ⁴ + c5*τ⁵
# 
# Method: use L(τ)/τ² = c2 + c3*τ + c4*τ² + c5*τ³
#         Extrapolate τ -> 0 using Richardson

def richardson_extrapolate(data, power):
    """
    Extrapolate L(τ)/τ^power to find c_{power}.
    Uses Neville's algorithm for Richardson extrapolation.
    """
    # Build sequence of L(τ)/τ^power
    mp.mp.dps = 150
    
    n = len(data)
    # Table for Richardson extrapolation
    R = [[mp.mpf('0')] * n for _ in range(n)]
    
    for i, (tau, L) in enumerate(data):
        tau_mp = mp.mpf(str(tau))
        L_mp = mp.mpf(str(L))
        R[i][0] = L_mp / (tau_mp**power)
    
    # Neville's algorithm
    for j in range(1, n):
        for i in range(n - j):
            tau_i = mp.mpf(str(data[i][0]))
            tau_j = mp.mpf(str(data[i+j][0]))
            R[i][j] = (tau_j**j * R[i][j-1] - tau_i**j * R[i+1][j-1]) / (tau_j**j - tau_i**j)
    
    return R[0][n-1]

# Extrapolate c2 from L(τ)/τ²
print("\nExtracting c2 from L(τ)/τ²:")
c2_rich = richardson_extrapolate(data, 2)
print(f"c2 (Richardson) = {mp.nstr(c2_rich, 30)}")
print(f"c2 (claimed)    = {mp.nstr(-mp.pi**2/9, 30)}")
print(f"Diff: {mp.nstr(c2_rich - (-mp.pi**2/9), 20)}")

# To get c3, compute (L/τ² - c2)/τ
print("\nExtracting c3 from (L(τ)/τ² - c2)/τ:")

def extract_c3(data, c2):
    mp.mp.dps = 150
    c2_mp = mp.mpf(str(c2))
    
    # Build sequence of (L/τ² - c2)/τ
    values = []
    for tau, L in data:
        tau_mp = mp.mpf(str(tau))
        L_mp = mp.mpf(str(L))
        val = (L_mp / tau_mp**2 - c2_mp) / tau_mp
        values.append((tau, val))
    
    # Extrapolate to τ -> 0
    n = len(values)
    R = [[mp.mpf('0')] * n for _ in range(n)]
    
    for i, (tau, val) in enumerate(values):
        R[i][0] = val
    
    for j in range(1, n):
        for i in range(n - j):
            tau_i = mp.mpf(str(data[i][0]))
            tau_j = mp.mpf(str(data[i+j][0]))
            R[i][j] = (tau_j**j * R[i][j-1] - tau_i**j * R[i+1][j-1]) / (tau_j**j - tau_i**j)
    
    return R[0][n-1]

c3_rich = extract_c3(data, c2_rich)
print(f"c3 (Richardson) = {mp.nstr(c3_rich, 30)}")
print(f"c3 (claimed)    = {mp.nstr(-16*mp.zeta(3)/9, 30)}")
print(f"Diff: {mp.nstr(c3_rich - (-16*mp.zeta(3)/9), 20)}")