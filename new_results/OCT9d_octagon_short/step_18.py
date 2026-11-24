import time
t0=time.time()
# Build [G|D] matrices and verify stored solutions reproduce TW3 exactly (PROVED check)
Gargs, Dargs = rhs_terms_235()
sols7, negw7, logprods7 = rhs_terms_2357()

def col_Li3(g, SS):
    return SS.coord(S_Li3(g))
def col_prod(w, p, SS):
    return SS.coord(S_mLi2logp(w, p))
def col_logprod(trip, SS):
    return SS.coord(S_logprod(trip))

# edges 1-3: columns 99 + 27
cols235 = []
for g in Gargs: cols235.append(col_Li3(g, SS235))
for (w,p) in Dargs: cols235.append(col_prod(w, p, SS235))
M235 = sp.Matrix.hstack(*cols235)
print("M235", M235.rows, M235.cols, "rank...", time.time()-t0)
rk = M235.rank(); print("rank =", rk, f"({time.time()-t0:.1f}s)")
for e in (1,2,3):
    sol = sp.Matrix([sp.Rational(s) for s in D['s123'][e]['sol']])
    diff = M235*sol - TW3[e]
    print(f"edge {e}: [G|D] sol - TW3 = 0 ?", diff == sp.zeros(27,1))
