
import mpmath as mp
mp.mp.dps = 130
import json
print("dps:", mp.mp.dps)
def Cm_formula(m):
    h = mp.pi/(2*m)
    s = mp.sin(h)
    return 48*m*s**4*mp.cos(h)*(3*s**2-1)
z3 = json.load(open('z3coeff_m2_8.json'))
print(f"{chr(39)}m{chr(39)}: formula vs json C_m  (residual)")
for m in range(2,9):
    f = Cm_formula(m)
    j = mp.mpf(z3[str(m)]['coeff'])
    print(f"m={m}: resid = {mp.nstr(abs(f-j),3)}   C_m={mp.nstr(f,30)}")
# Known closed forms
print()
checks = {2:6*mp.sqrt(2), 3:-9*mp.sqrt(3)/8, 6:(189*mp.sqrt(6)-333*mp.sqrt(2))/8}
for m,v in checks.items():
    print(f"m={m}: formula - known = {mp.nstr(Cm_formula(m)-v,3)}")
# m=8 form: -9g^7+96g^5-336g^3+384g, g=2cos(pi/16)
g = 2*mp.cos(mp.pi/16)
c8 = -9*g**7+96*g**5-336*g**3+384*g
print("m=8: formula - poly(g) =", mp.nstr(Cm_formula(8)-c8,3))
# Predict m=9,10,12
print()
print("PREDICTIONS:")
for m in [9,10,12]:
    print(f"  C_{m} = {mp.nstr(Cm_formula(m),120)}")
    with open(f'pred_m{m}.txt','w') as fo: fo.write(mp.nstr(Cm_formula(m),120))
