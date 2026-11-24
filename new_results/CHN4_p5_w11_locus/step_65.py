import pickle, mpmath as mp
up=pickle.load(open('gamma_up.pkl','rb'))
down=pickle.load(open('gamma_down.pkl','rb'))
base=pickle.load(open('gamma_pt1.pkl','rb'))
# assemble Gamma points: (a1,a2,a4,a5)
pts=[]
for r in down: pts.append(r)              # a1=0.0005..0.005
pts.append((base['a1'],base['a2'],base['a4'],base['a5']))  # 0.01
for r in up: pts.append(r)                # 0.02,0.05,0.1,0.2
pts=sorted(pts, key=lambda r: float(r[0]))
print("Gamma points for W12 batch:")
sel=[pts[i] for i in [0,2,4,6,7,8]]  # spread: 0.0005,0.002,0.01,0.05,0.1,0.2
for r in pts: print("  a1=",mp.nstr(r[0],6)," a2=",mp.nstr(r[1],6)," a4=",mp.nstr(r[2],6)," a5=",mp.nstr(r[3],6))
print("\nselected for W12:")
for r in sel: print("  a1=",mp.nstr(r[0],6))
pickle.dump(sel, open('gamma_batch.pkl','wb'))
print("saved gamma_batch.pkl with", len(sel), "points")
