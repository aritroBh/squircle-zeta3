import sympy as sp, pickle, os, math
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
r5 = sp.sqrt(5); rho = sp.sqrt(10-2*r5)
phi = (1+r5)/2
t = sp.Symbol('t')
h  = (r5-1)/4; y0 = (r5-2)/rho
# exact edge identities
def pq(e):
    e = sp.expand(e)
    p = sp.simplify((e + e.subs(r5,-r5))/2)
    q = sp.simplify((e - e.subs(r5,-r5))/(2*r5))
    return sp.Rational(p), sp.Rational(q)
assert sp.expand(2*phi*h - 1) == 0                       # h = 1/(2 phi)
L2 = sp.expand(4*y0**2)                                   # L^2
assert sp.simplify(L2 - (sp.Rational(25,10) - sp.Rational(11,10)*r5)) == 0
print("L^2 = 4*y0^2 = (25-11 sqrt5)/10  [EXACT]")
h2oL2 = sp.simplify(h**2/L2)
assert sp.simplify(h2oL2 - (sp.Rational(5,4)+r5/2)) == 0
print("h^2/L^2 = (5+2 sqrt5)/4  [EXACT]")
n = [ ( (1+r5)/4, rho/4 ), ( (r5-1)/4, rho*(1+r5)/8 ), ( (1-r5)/4, rho*(1+r5)/8 ),
      ( -(1+r5)/4, rho/4 ), ( sp.Integer(-1), sp.Integer(0) ) ]
V0 = (-h, y0); V1 = (-h, -y0)
y  = (V0[0] + t*(V1[0]-V0[0]), V0[1] + t*(V1[1]-V0[1]))
yp = (sp.diff(y[0],t), sp.diff(y[1],t))
yy  = sp.expand(y[0]**2+y[1]**2); yyp = sp.expand(y[0]*yp[0]+y[1]*yp[1]); yp2 = sp.expand(yp[0]**2+yp[1]**2)
assert sp.simplify(yp2 - L2) == 0                                  # |y'| = L
assert sp.simplify(yyp/L2 + (1-2*t)/2) == 0                        # y.y'/L^2 = -(1-2t)/2
assert sp.simplify(yy/L2 - (h2oL2 + (t-sp.Rational(1,2))**2)) == 0 # |y|^2/L^2
print("|y'|^2=L^2, y.y'/L^2=-(1-2t)/2, |y|^2/L^2=h^2/L^2+(t-1/2)^2  [EXACT]")
# affine term values, DEC1 order [a2, a3, a1, a4, a5]
sig_idx = [(1,-1),(2,1),(0,-1),(3,1),(4,1)]  # (my_i_index, sigma)  from previous step
avals = []
for i in range(5):
    ex = y[0]*n[i][0] + y[1]*n[i][1]
    avals.append(sp.expand(ex))
apol = []
appol = []
for mi, sg in sig_idx:
    a = sp.expand(sg*avals[mi])
    a0 = sp.expand(a.subs(t,0)); a1 = sp.expand(a.subs(t,1)-a0)
    apol.append([a0,a1]); appol.append([a1, sp.Integer(0)])
# cross-check with DEC1's constants
A0 = 1/(2*phi**2); A1 = A0; B1 = A0; B4 = 1/(2*phi); B5 = B4
kap1 = phi-1; kap4 = 1/phi-1
dec1 = [[sp.Integer(0),A0],[A1,-A1],[B1,B1*kap1],[B4,B4*kap4],[B5,sp.Integer(0)]]
for i in range(5):
    for j in range(2):
        assert sp.simplify(apol[i][j]-sp.expand(dec1[i][j]))==0, (i,j)
print("term values match DEC1 exactly: a2=A0 t, a3=A1(1-t), a1=B1(1+t/phi), a4=B4(1-t/phi^2), a5=B5")
# sum = 1 exactly
assert sp.expand(sum(sp.expand(a[0]+a[1]*t) for a in apol) - 1) == 0
print("sum a_i(t) = 1  [EXACT]  (edge lies in dK1)")
# derivatives
for i in range(5):
    assert sp.expand(appol[i][0] - sp.diff(apol[i][0]+apol[i][1]*t, t)) == 0
# geometry identity for P formula: rho^2=10-2r5 used implicitly; verify y0^2 = h^2*cot^2? skip
pickle.dump({'apol':apol,'appol':appol,'h':h,'y0':y0,'L2':L2,'h2oL2':h2oL2}, open('geom_exact.pkl','wb'))
with open('progress.log','a') as f:
    f.write("step_01/02: exact geometry derived (normals, V0,V1, term values, h, L^2=(25-11sqrt5)/10, h2oL2); all identities EXACT; matches DEC1 geometry.pkl\n")
print("saved geom_exact.pkl")
