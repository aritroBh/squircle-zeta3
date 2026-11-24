"""step_02.py -- identify K_1 = polar of zonotope Z = sum_j [-n_j,n_j] for models (a),(b)."""
import os, itertools, json
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 60
import numpy as np
from scipy.spatial import ConvexHull

sq2 = mp.sqrt(2); sq3 = mp.sqrt(3)
N_A = [[1/sq2,1/sq2,0],[1/sq2,-1/sq2,0],[1/sq2,0,1/sq2],[1/sq2,0,-1/sq2],[0,1/sq2,1/sq2],[0,1/sq2,-1/sq2]]
N_B = [[1,0,0],[0,1,0],[0,0,1],[1/sq3,1/sq3,1/sq3],[1/sq3,1/sq3,-1/sq3],[1/sq3,-1/sq3,1/sq3],[-1/sq3,1/sq3,1/sq3]]

def zonotope(normals):
    m = len(normals)
    pts = []
    for sig in itertools.product([1,-1], repeat=m):
        p = [sum(sig[j]*normals[j][i] for j in range(m)) for i in range(3)]
        pts.append(p)
    pts = np.array(pts, float)
    hull = ConvexHull(pts)
    # unique extreme vertices
    V = np.unique(np.round(pts[hull.vertices], 12), axis=0)
    # group triangular facets by plane -> polygonal faces
    faces = {}
    for simp, nrm in zip(hull.simplices, hull.equations):
        key = tuple(np.round(nrm, 9))
        faces.setdefault(key, set()).update(simp)
    return pts, hull, V, faces

def polar_faces(normals, Zfaces):
    """K_1 facets: plane <y, vZ> = 1 for each vertex vZ of Z (extreme)."""
    pass

def analyze(name, normals):
    print("="*70); print(name); print("="*70)
    pts, hull, V, faces = zonotope(normals)
    print(f"Z: {len(V)} extreme vertices, {len(faces)} polygonal facets")
    # face size histogram
    from collections import Counter
    hist = Counter(len(v) for v in faces.values())
    print("Z facet size histogram:", dict(hist))
    # K1 vertices = polar of Z facets: y = n/(-d) for plane n.x + d = 0  (d<0 since origin inside)
    K1V = []
    for key, verts in faces.items():
        eq = np.array(key)  # n.x + d = 0 with |n|=1
        n, d = eq[:3], eq[3]
        y = -n/d   # polar point: <y,x> <= 1 with equality on plane
        K1V.append(y)
    K1V = np.array(K1V)
    print(f"K1: {len(K1V)} vertices")
    # K1 facets: for each Z-vertex v (extreme), plane <y,v>=1; find K1 vertices on it
    tol = 1e-9
    facets = []
    for v in V:
        on = [i for i,y in enumerate(K1V) if abs(np.dot(y,v)-1) < tol]
        facets.append(sorted(on))
    fhist = Counter(len(f) for f in facets)
    print("K1 facet size histogram:", dict(fhist))
    # edges: pairs of K1 vertices sharing >=2 facets
    from collections import defaultdict
    vfacets = defaultdict(list)
    for fi,f in enumerate(facets):
        for vi in f: vfacets[vi].append(fi)
    edges = set()
    for f in facets:
        for a,b in itertools.combinations(f,2):
            shared = set(vfacets[a])&set(vfacets[b])
            if len(shared) >= 2: edges.add((min(a,b),max(a,b)))
    edges = sorted(edges)
    print(f"K1: {len(edges)} edges, Euler check V-E+F = {len(K1V)-len(edges)+len(facets)}")
    # surface area
    def tri_area(a,b,c):
        return 0.5*np.linalg.norm(np.cross(b-a, c-a))
    S = 0.0
    for f in facets:
        # fan triangulate (vertices may be unordered; order them around centroid)
        ptsf = K1V[f]
        c = ptsf.mean(0)
        nrm = V[[i for i,v in enumerate(V) if False]] # placeholder
        # normal from first 3 noncollinear
        nn = np.cross(ptsf[1]-ptsf[0], ptsf[2]-ptsf[0]); nn/=np.linalg.norm(nn)
        ax1 = ptsf[0]-c; ax1/=np.linalg.norm(ax1); ax2 = np.cross(nn, ax1)
        angs = np.arctan2((ptsf-c)@ax2, (ptsf-c)@ax1)
        order = np.argsort(angs)
        for i in range(1,len(f)-1):
            S += tri_area(ptsf[order[0]], ptsf[order[i]], ptsf[order[i+1]])
    print(f"S(1) = {mp.nstr(mpf(float(S)),40)}")
    # edge data: lengths, dihedrals, adjacent face sizes
    # face outward unit normals of K1: for facet from Z-vertex v: normal = v/|v|
    fnrm = []
    for i, f in enumerate(facets):
        v = None
        # recover the Z vertex used
        fnrm.append(None)
    # rebuild: facets[i] came from V[i'] in same loop order; redo mapping
    facets2 = []
    for v in V:
        on = [i for i,y in enumerate(K1V) if abs(np.dot(y,v)-1) < tol]
        facets2.append((sorted(on), v/np.linalg.norm(v)))
    vfacets2 = defaultdict(list)
    for fi,(f,nr) in enumerate(facets2):
        for vi in f: vfacets2[vi].append(fi)
    edata = []
    for (a,b) in edges:
        shared = sorted(set(vfacets2[a])&set(vfacets2[b]))
        f1, f2 = shared[0], shared[1]
        n1, n2 = facets2[f1][1], facets2[f2][1]
        dih = np.pi - np.arccos(np.clip(np.dot(n1,n2),-1,1))
        L = np.linalg.norm(K1V[a]-K1V[b])
        edata.append((a,b,float(L),float(dih),len(facets2[f1][0]),len(facets2[f2][0])))
    # group edges by (L, dihedral, face sizes)
    groups = defaultdict(list)
    for e in edata:
        key = (round(e[2],9), round(e[3],9), tuple(sorted((e[4],e[5]))))
        groups[key].append(e)
    print("edge orbits (L, dihedral, faces):")
    for k,v in sorted(groups.items()):
        print(f"   count={len(v):3d}  L={k[0]:.12f}  dih={k[1]:.12f}  faces={k[2]}")
    return K1V, facets2, edata, float(S)

KA = analyze("MODEL (a): 6 face-diagonal normals", N_A)
KB = analyze("MODEL (b): 3 axes + 4 body diagonals", N_B)
np.save("K1V_A.npy", KA[0]); np.save("K1V_B.npy", KB[0])
