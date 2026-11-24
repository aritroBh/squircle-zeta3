# Planner check of K1g: P^(5)(1) of the golden hexagon V=(1,0),(phi,1),(0,1), w=(2,3,5), by my own polar code + 5th-difference + Richardson.
import mpmath as mp
mp.mp.dps = 70
phi = (1+mp.sqrt(5))/2
V = [(mp.mpf(1),mp.mpf(0)), (phi,mp.mpf(1)), (mp.mpf(0),mp.mpf(1))]; w = [2,3,5]
cr = lambda P,Q: P[0]*Q[1]-P[1]*Q[0]
A = mp.matrix([[abs(cr(V[j],V[q])) for j in range(3)] for q in range(3)])  # N(V_q)=1 -> w_j k_j
c = mp.lu_solve(A, mp.matrix([1,1,1])); k = [c[j]/w[j] for j in range(3)]
def P(p):
    cc = [w[j]*k[j]**p for j in range(3)]
    def g(th):
        s, co = mp.sin(th), mp.cos(th); F = 0; Fp = 0
        for j in range(3):
            sj = V[j][0]*s - V[j][1]*co; dsj = V[j][0]*co + V[j][1]*s
            F += cc[j]*abs(sj)**p; Fp += cc[j]*p*abs(sj)**(p-1)*mp.sign(sj)*dsj
        rho = F**(-1/p); drho = -(1/p)*F**(-1/p-1)*Fp
        return mp.sqrt(rho**2+drho**2)
    br = [mp.atan2(v[1],v[0]) for v in V] + [mp.pi]
    return 2*sum(mp.quad(g,[br[i],br[i+1]]) for i in range(3))
def d5(h):
    f = {j: P(1+j*h) for j in (-3,-2,-1,1,2,3)}
    return (-f[-3]+4*f[-2]-5*f[-1]+5*f[1]-4*f[2]+f[3])/(2*h**5)
hs = [mp.mpf('4e-4'), mp.mpf('2e-4'), mp.mpf('1e-4')]
D = [d5(h) for h in hs]
R1 = [(4*D[i+1]-D[i])/3 for i in range(2)]; R2 = (16*R1[1]-R1[0])/15
print("d5:", [mp.nstr(x,18) for x in D]); print("Richardson:", mp.nstr(R2,20))
print("K1g formula: -78.3585998739646021275  diff:", mp.nstr(R2-mp.mpf('-78.3585998739646021275'),3))
