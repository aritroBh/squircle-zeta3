import sys, json, importlib.util
import mpmath as mp

# --- load mzv.py from E5c ---
spec = importlib.util.spec_from_file_location("mzv", "../E5c_nonlocal_z5/mzv.py")
mzv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mzv)

mp.mp.dps = 360
d = json.load(open('../E5d_nonlocal_pslq/b0_p5_convB_350.json'))
e1 = mp.mpf(d['e1']['P5']); e2 = mp.mpf(d['e2']['P5'])
print("e1 - e2 =", mp.nstr(e1-e2, 10))
print("rel diff =", mp.nstr(abs(e1-e2)/abs(e1), 5))

B = mzv.basis(360)
print("basis size:", len(B))
assert len(B) == 20
# after basis() call, dps may be reset; ensure
mp.mp.dps = 360
sq5 = mp.sqrt(5)
X = e1/sq5
print("X =", mp.nstr(X, 40))

# sanity: compare with E5d 187-digit value
old = "-852.255256336156211501665304739892950203528360898306405873984692613876904867441"
print("prefix match with E5d 87-digit fragment:", d['e1']['P5'].startswith(old))

vals = {"X": mp.nstr(X, 350)}
for (lab, v, w) in B:
    vals[lab] = mp.nstr(v, 350)
json.dump({"labels": [b[0] for b in B], "weights": [b[2] for b in B], "vals": vals},
          open("consts360.json","w"), indent=1)
print("saved consts360.json")
for i,(lab,v,w) in enumerate(B):
    print(f"{i+1:3d} w{w} {lab:16s} {mp.nstr(v,30)}")
print("  X            ", mp.nstr(X,30))
