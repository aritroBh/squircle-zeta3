"""Convention dependence: same K_1 as the parallelogram (1,1/2),(-1/3,1), weights w=(2,3) held fixed:
K_p = {sum_j w_j |<v_j/w_j, x>|^p <= 1}."""
import mpmath as mp
from jet import perimeter_jet
from geom import edges, edge_data
mp.mp.dps = 60
R = lambda a, b=1: mp.mpf(a)/b
V = [(R(1),R(1,2)),(R(-1,3),R(1))]; w = [R(2), R(3)]
U = [(V[j][0]/w[j], V[j][1]/w[j]) for j in range(2)]
P = perimeter_jet(U, W=w); P0 = perimeter_jet(V)
E = mp.fsum(edge_data(V, *e)[1]**2/edge_data(V, *e)[0] for e in edges(V))
print("weighted:   P', P'', P''' =", [mp.nstr(x, 25) for x in (P[1], 2*P[2], 6*P[3])])
print("unweighted: P', P'', P''' =", [mp.nstr(x, 25) for x in (P0[1], 2*P0[2], 6*P0[3])])
print("printed law: 0,", mp.nstr(mp.pi**2*E/3, 25), mp.nstr(6*E*(mp.zeta(3)-mp.pi**2/3), 25))
