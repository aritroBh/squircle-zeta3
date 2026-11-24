
import sympy as sp

# Symbols: vanishing sines/cosines and their logs
s0, s1, c0, c1 = sp.symbols('s0 s1 c0 c1')
L0, L1 = sp.symbols('L0 L1')          # the two vanishing logs (l0, l1)
# Reduced sums over NON-vanishing indices N
tau1, tau2, tau3 = sp.symbols('tau1 tau2 tau3')      # sum_{N} s_i l_i^k
gam1, gam2, gam3 = sp.symbols('gam1 gam2 gam3')      # sum_{N} c_i(l_i^k + k l_i^{k-1})
S0 = sp.Symbol('S0')        # full sum of sines  = s0+s1+sum_N s_i
Sc = sp.Symbol('Sc')        # full sum of cosines = c0+c1+sum_N c_i
lS0 = sp.Symbol('lS0')      # log(S0)

# Full T_k (k=0..3) and Tp_k as polynomials in L0,L1
T0 = S0
T1 = tau1 + s0*L0 + s1*L1
T2 = tau2 + s0*L0**2 + s1*L1**2
T3 = tau3 + s0*L0**3 + s1*L1**3
Tp0 = Sc
Tp1 = gam1 + c0*(L0 + 1) + c1*(L1 + 1)
Tp2 = gam2 + c0*(L0**2 + 2*L0) + c1*(L1**2 + 2*L1)
Tp3 = gam3 + c0*(L0**3 + 3*L0**2) + c1*(L1**3 + 3*L1**2)

T = [T0,T1,T2,T3]; Tp=[Tp0,Tp1,Tp2,Tp3]
S0p = Tp0
U1 = T[1]/S0
U2 = T[2]/(2*S0) - T[1]**2/(2*S0**2)
U3 = T[3]/(6*S0) - T[1]*T[2]/(2*S0**2) + T[1]**3/(3*S0**3)
U1p = Tp[1]/S0 - T[1]*S0p/S0**2
U2p = Tp[2]/(2*S0) - T[2]*S0p/(2*S0**2) - T[1]*Tp[1]/S0**2 + T[1]**2*S0p/S0**3
U3p = (Tp[3]/(6*S0) - T[3]*S0p/(6*S0**2)
       - (Tp[1]*T[2]/(2*S0**2) + T[1]*Tp[2]/(2*S0**2) - T[1]*T[2]*S0p/S0**3)
       + (T[1]**2*Tp[1]/S0**3 - T[1]**3*S0p/S0**4))
L = S0p/S0
R1 = lS0 - U1; R2 = -lS0 + U1 - U2; R3 = lS0 - U1 + U2 - U3
Q0 = -L
Q1 = L - U1p; Q2 = -L + U1p - U2p; Q3 = L - U1p + U2p - U3p
rho1 = R1; rho2 = R2 + R1**2/2; rho3 = R3 + R1*R2 + R1**3/6
W1 = 2*Q0*Q1; W2 = 2*Q0*Q2 + Q1**2; W3 = 2*Q0*Q3 + 2*Q1*Q2
B = 1 + Q0**2
E1 = W1/(2*B); E2 = W2/(2*B) - W1**2/(8*B**2); E3 = W3/(2*B) - W1*W2/(4*B**2) + W1**3/(16*B**3)
D0 = (1/S0)*sp.sqrt(B)
D3 = D0*(rho3 + rho2*E1 + rho1*E2 + E3)
D3 = sp.expand(D3)
print("D3 expanded. Number of terms:", len(D3.as_ordered_terms()))
# Extract coefficient of L0^2 L1 and L0 L1^2
G21 = sp.expand(D3).coeff(L0,2).coeff(L1,1)
G12 = sp.expand(D3).coeff(L0,1).coeff(L1,2)
print("G21 extracted, terms:", len(sp.expand(G21).as_ordered_terms()))
print("G12 extracted, terms:", len(sp.expand(G12).as_ordered_terms()))
