
import pickle, sympy as sp, glob
# verify Q8, Q9 transcription
Qpolys = pickle.load(open('Qpolys.pkl','rb'))
u = sp.symbols('u')
z = {i: sp.Symbol(f'z{i}') for i in range(2,10)}
Q8_txt = "-(2891923/12288 + 495965*u/256 + 60889*u**2/48 + 127*u**3)*z8 + (141/16 + 50*u + 16*u**2)*z3*z5 + (169/32 + 27*u + 8*u**2)*z2*z3**2"
Q9_txt = ("-(22805/192 + 3907*u/4 + 1916*u**2/3 + 64*u**3)*z9 - (8343/64 + 4163*u/4 + 660*u**2 + 64*u**3)*z2*z7 "
          "- (60093/640 + 29741*u/40 + 2334*u**2/5 + 224*u**3/5)*z2**2*z5 - (263713/4480 + 130973*u/280 + 2062*u**2/7 + 992*u**3/35)*z2**3*z3 "
          "+ (169/96 + 9*u + 8*u**2/3)*z3**3")
loc = {f'z{i}': z[i] for i in range(2,10)} | {'u': u}
Q8 = sp.sympify(Q8_txt.replace('z8','z8').replace('z3*z5','z3*z5'), locals=loc)
Q9 = sp.sympify(Q9_txt, locals=loc)
t8 = sp.sympify(Qpolys[8], locals=loc)
t9 = sp.sympify(Qpolys[9], locals=loc)
print("Q8 match:", sp.expand(Q8 - t8) == 0)
print("Q9 match:", sp.expand(Q9 - t9) == 0)

# count I certificates
icerts = sorted(glob.glob('cert_I_*.pkl'))
print("I certificates:", len(icerts), [f.split('cert_I_')[1].replace('.pkl','') for f in icerts])
