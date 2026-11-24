src = '''"""obstruction32.py -- the antisymmetrised (3,2) cobracket of a weight-5 symbol,
projected to Lie_3(V) (x) Lie_2(V), for V of any dimension r >= 2 and any exact
coefficient ring (duck-typed: coefficients need +, *, ==; zero = the int/Fraction 0
or a ring zero).

Symbol: dict {word: coeff}, word = tuple of slot indices in {0..r-1}.
Psi(S) = sum_w s(w) * (w[:3] (x) w[3:]  -  w[2:] (x) w[:2])   in V^t3 (x) V^t2.
Pairing with A (tensor elt of V^t3) and B (tensor elt of V^t2):
   <A(x)B, Psi> = sum_{(u,v)} Psi[(u,v)] A[u] B[v].
Lie bases (Hall): Lie_2: [X_i,X_j], i<j.  Lie_3: [X_i,[X_j,X_k]], j<k, j<=i.
The obstruction vanishes iff all basis pairings vanish (PROVED symbol-level test).
Controls: symbols of Li_5 and of shuffle products pair to zero.
"""

def lie2_basis(r):
    """list of (name, {word:int}) for [X_i, X_j], i<j."""
    out = []
    for i in range(r):
        for j in range(i+1, r):
            out.append(((i, j), {(i, j): 1, (j, i): -1}))
    return out

def lie3_basis(r):
    """list of (name, {word:int}) for [X_i, [X_j, X_k]], j<k, j<=i."""
    out = []
    for j in range(r):
        for k in range(j+1, r):
            for i in range(j, r):
                out.append(((i, j, k), {(i, j, k): 1, (i, k, j): -1,
                                        (j, k, i): -1, (k, j, i): 1}))
    return out

def psi_tensor(S):
    """Psi as dict {(u3, v2): coeff} (integer-coefficient arithmetic on the coeff ring)."""
    Psi = {}
    for w, c in S.items():
        if len(w) != 5 or c == 0: continue
        k1 = (w[:3], w[3:])
        Psi[k1] = Psi.get(k1, 0) + c
        k2 = (w[2:], w[:2])
        Psi[k2] = Psi.get(k2, 0) - c
    return {k: v for k, v in Psi.items() if v != 0}

def pair(Psi, A, B):
    """<A (x) B, Psi>; A, B dicts word->int."""
    tot = 0
    for (u, v), c in Psi.items():
        av = A.get(u)
        if av is None: continue
        bv = B.get(v)
        if bv is None: continue
        tot = tot + c*av*bv
    return tot

def obstruction(S, r):
    """full table of pairings against Lie_3 (x) Lie_2 Hall bases.
    Returns (table, names): table[i][j] = pairing with (Lie3_i (x) Lie2_j)."""
    Psi = psi_tensor(S)
    L3 = lie3_basis(r); L2 = lie2_basis(r)
    table = [[pair(Psi, A, B) for (_, B) in L2] for (_, A) in L3]
    return table, ([n for (n, _) in L3], [n for (n, _) in L2])

def is_zero(table):
    return all(v == 0 for row in table for v in row)
'''
with open("obstruction32.py","w") as f: f.write(src)
print("written obstruction32.py")
