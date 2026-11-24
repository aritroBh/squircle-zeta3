from mpmath import mp, mpf, pi, sin, cos, tan, log, quad, polylog
mp.dps = 60

def Cl2(theta):
    return (polylog(2, mp.e**(1j*theta))).imag

def L1_closed(phi):
    # (1/2)[tan ln sin + ∫ ln cot] - phi/2 ; ∫_0^phi ln cot u du = (Cl2(pi-2phi)+Cl2(2phi))/2
    return mpf('0.5')*(tan(phi)*log(sin(phi)) + (Cl2(pi-2*phi)+Cl2(2*phi))/2) - phi/2

def L1_quad(phi):
    return quad(H1u, [0, phi])

for phi in [0.3, 0.8, 1.2, 1.5]:
    a = L1_quad(phi); b = L1_closed(phi)
    print(phi, mp.nstr(a-b, 5))
