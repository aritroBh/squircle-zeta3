# Referee items 1 + 3: symbolic checks (own code).
# Truncated power series in t are handled by hand (polynomial truncation) for speed.
import sympy as sp
t,s,H,S,T,tau = sp.symbols('t s H S T tau', real=True)
N = 4
def trunc(e, n=N):  # keep t^0..t^(n-1)
    e = sp.expand(e)
    return sum(e.coeff(t,k)*t**k for k in range(n))
def exp_series(x, n=N):  # exp(x) with x = O(t), truncated at t^n
    out, term = sp.Integer(1), sp.Integer(1)
    for k in range(1,n):
        term = trunc(term*x, n)/k
        out += term
    return trunc(out, n)
# --- B_n
a = sp.symbols('a1:5')
ser = exp_series(sum(t**k*a[k-1] for k in range(1,N+1)), N+1)
B = [sp.expand(ser.coeff(t,n)) for n in range(N+1)]
print("B_n:", B)
print("B2 check:", sp.simplify(B[2]-(a[1]+a[0]**2/2)), " B3 check:", sp.simplify(B[3]-(a[2]+a[0]*a[1]+a[0]**3/6)))
ad = {(k,i): sp.Symbol(f'a{k}_d{i}') for k in range(1,N) for i in range(0,N-k)}
g = sp.symbols('g1:4')
for sg in (+1,-1):
    absS = sg*s
    cosser = sum((-1)**j*(t*s)**(2*j)/sp.factorial(2*j) for j in range(3))
    sinser = sum((-1)**j*(t*s)**(2*j+1)/sp.factorial(2*j+1) for j in range(3))
    Mser = trunc(H*cosser + S*sinser + T*sg*sinser)    # |sin(ts)| = sg*sin(ts) for small t>0
    Pser = exp_series(sum(t**k*sum(ad[(k,i)]*(t*s)**i/sp.factorial(i) for i in range(0,N-k)) for k in range(1,N)))
    Gser = exp_series(sum(t**k*g[k-1] for k in range(1,N))) - 1
    prod = trunc(Mser*Pser*Gser)
    Phi = [sp.expand(prod.coeff(t,n)) for n in range(N)]
    M1 = S*s+T*absS; M2 = -H*s**2/2
    P1 = ad[(1,0)]; P2 = ad[(2,0)]+ad[(1,0)]**2/2+ad[(1,1)]*s
    G1=g[0]; G2=g[1]+g[0]**2/2; G3=g[2]+g[0]*g[1]+g[0]**3/6
    claim = {1: H*G1, 2: H*G2+(M1+H*P1)*G1, 3: H*G3+(M1+H*P1)*G2+(M2+M1*P1+H*P2)*G1}
    print(f"s sign {sg:+d}: Phi_0 = {Phi[0]}")
    for n in (1,2,3):
        print(f"   Phi_{n} - claimed = {sp.expand(Phi[n]-claim[n])}")
    print("   M_a:", [sp.factor(Mser.coeff(t,k)) for k in range(N)], "  P_2:", sp.expand(Pser.coeff(t,2)))
    print("   P_3 (not used in c4, for reference):", sp.expand(Pser.coeff(t,3)))

# --- Proposition 4: exact derivative algebra. c=cos, sn=sin, E2,E3,E4 symbols with
#     D c = -sn, D sn = c, D E_i = (psi_i'/tau) E_i.
c, sn, E2, E3, E4 = sp.symbols('c sn E2 E3 E4', positive=True)
M = c+sn
psi = {E2: -2*sn, E3: -2*c, E4: -2*M}
def D(e):
    e = sp.sympify(e)
    r = sp.diff(e,c)*(-sn) + sp.diff(e,sn)*c
    for Ei,p in psi.items():
        r += sp.diff(e,Ei)*D0(p)/tau*Ei
    return r
def D0(e):  # derivative of a polynomial in c, sn only
    return sp.diff(e,c)*(-sn)+sp.diff(e,sn)*c
eps = E2+E3+E4
ellp = D(eps)/(1+eps)
ellpp = D(ellp)
Mp = D0(M)
Q = 2*Mp*ellp + M*ellpp + tau*M*ellp**2
lhs = tau**2*(1+eps)**2*Q
rhs = (1+eps)*sum(Ei*(tau*(2*Mp*D0(p)-M*p)+M*D0(p)**2) for Ei,p in psi.items()) \
      - (1-tau)*M*(sum(D0(p)*Ei for Ei,p in psi.items()))**2
diffe = sp.simplify(sp.together(lhs-rhs))
print("Prop4 identity tau^2(1+eps)^2 Q - RHS =", diffe)
# Phi + Phi'' = tau f Q with f = exp(tau*ell): Phi = M f, f' = tau ell' f, f'' = (tau ell'' + tau^2 ell'^2) f
fp = tau*ellp; fpp = tau*ellpp + tau**2*ellp**2       # divided by f
PhiPhi = M + (D0(D0(M)) + 2*Mp*fp + M*fpp)            # (Phi+Phi'')/f
print("(Phi+Phi'')/f - tau Q =", sp.simplify(sp.together(PhiPhi - tau*Q)))
for name,p in (('i=2',psi[E2]),('i=3',psi[E3]),('i=4',psi[E4])):
    print(name, " psi' =", D0(p), "  2M'psi'-M psi =", sp.expand(2*Mp*D0(p)-M*p))
# CS step: (1+eps)X - (1-tau)M Y^2 >= sum E_i[tau(1+eps)A_i + M psi_i'^2 (1+tau eps)]  given Y^2 <= eps*sum psi'^2 E
Ai = {Ei: 2*Mp*D0(p)-M*p for Ei,p in psi.items()}
lower = sum(Ei*(tau*(1+eps)*Ai[Ei] + M*D0(p)**2*(1+tau*eps)) for Ei,p in psi.items())
gap = sp.expand(rhs - lower)
target = sp.expand((1-tau)*M*(eps*sum(D0(p)**2*Ei for Ei,p in psi.items()) - (sum(D0(p)*Ei for Ei,p in psi.items()))**2))
print("RHS - lower bound - (1-tau)M[eps*sum psi'^2E - (sum psi'E)^2] =", sp.simplify(gap-target))
