import pickle, sys
sys.path.insert(0,'.')
import ftfield as F
g = pickle.load(open('geom14.pkl','rb'))
L = pickle.load(open('letters14.pkl','rb'))
print("geom14 keys:", list(g.keys()))
print("h =", g['h'], " eval", F.feval(g['h'],30))
print("L^2 =", g['L'])
print("A0 =", g['A0'], " A1 =", g['A1'])
print("B =", {k:g['B'][k] for k in g['B']})
print("X =", {k:g['X'][k] for k in g['X']})
for i in [1,2,5,6]:
    print(f"  x{i} eval {F.feval(g['X'][i],25)}   1-x{i} eval {F.feval(F.fsub(F.ONE,g['X'][i]),25)}  B{i} eval {F.feval(g['B'][i],25)}")
print()
print("letters14 keys:", list(L.keys()))
for k in L:
    v = L[k]
    s = str(v)
    print(" ", k, "=", s[:400])
