import numpy as np, mpmath as mp
mp.mp.dps = 50
WORK='research/08_theory/K4d_identify'
VA = np.load(WORK+'/K1V_A.npy')
print("K1V_A shape:", VA.shape)
s2 = mp.sqrt(2)
NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
# Identify vertices exactly
print("\nVertices (float) and candidate exact forms:")
for i,v in enumerate(VA):
    print(i, np.round(v,8), " |v|=%.6f"%np.linalg.norm(v))
# term values at each vertex: a_j = <n_j, V>
print("\nTerm values a_j(V) for each vertex (6 normals):")
for i,v in enumerate(VA):
    av = [float(np.dot(np.array([float(x) for x in n]), v)) for n in NA]
    print(i, np.round(av,6), " sum|.|=%.6f"%sum(abs(x) for x in av))
