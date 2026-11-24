
import pickle, sys
import sympy as sp
import mpmath as mp
mp.mp.dps = 50
base = "<project>/research/08_theory/"
sys.path.insert(0, base + "E5c_nonlocal_z5")
import symeng as SE

t = sp.Symbol('t')
el0, el1, W1, W4, L2, LF = sp.symbols('el0 el1 W1 W4 L2 LF')
E5d = sp.sympify(pickle.load(open("E5_dec_subs.pkl","rb")),
                 locals={'t':t,'el0':el0,'el1':el1,'W1':W1,'W4':W4,'L2':L2,'LF':LF})
E5f = sp.lambdify([t, el0, el1, W1, W4, L2, LF], E5d, modules='mpmath')

sqrt5 = mp.sqrt(5); phi = (1+mp.sqrt(5))/2
h = 1/(2*phi); R = 1/(2*phi**2*mp.sin(mp.pi/5)); Ry = R/(2*phi)
V0 = [-h, Ry]; V1 = [-h, -Ry]
M = [[mp.cos(j*mp.pi/5), mp.sin(j*mp.pi/5)] for j in range(5)]
P = SE.edge_params(M, V0, V1)
L = mp.sqrt((25-11*sqrt5)/10)

maxd = mp.mpf(0)
for tv in ['0.13','0.3','0.5','0.71','0.88']:
    tv = mp.mpf(tv)
    val_ref = SE.E5_sym(tv, M, P, conv='B')
    val_mine = L*E5f(tv, mp.log(tv), mp.log(1-tv), mp.log(1+tv/phi), mp.log(1-tv/phi**2), mp.log(2), mp.log(phi))
    d = abs(val_ref - val_mine)
    maxd = max(maxd, d)
    print(f"t={tv}: symeng={mp.nstr(val_ref,18)}  mine={mp.nstr(val_mine,18)}  diff={mp.nstr(d,3)}")
print("\nMAX DIFF:", mp.nstr(maxd, 4))
# normalization check: P5 = 1200 * int E5
I = mp.quad(lambda x: SE.E5_sym(x, M, P, conv='B'), [0,1])
P5gold = mp.mpf(pickle.load(open("p5_gold.pkl","rb"))['P5'])
print("\nint E5 =", mp.nstr(I, 30))
print("1200*I =", mp.nstr(1200*I, 30))
print("P5gold =", mp.nstr(P5gold, 30))
print("diff:", mp.nstr(abs(1200*I - P5gold), 3))
