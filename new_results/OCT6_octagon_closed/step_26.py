import pickle, mpmath as mp, sympy as sp
import mixed_engine as ME
mp.mp.dps = 730
num = {sp.symbols(k): mp.mpf(v) for k,v in pickle.load(open('num745.pkl','rb')).items()}
S = pickle.load(open('S_assembly.pkl','rb'))
r2m=mp.sqrt(2); um=mp.sqrt(2+r2m)
def ev(e): return mp.mpf(sp.N(e.subs(num, simultaneous=True), 720))
t0=__import__('time').time()
v_mix = ev(S['S_mix']); v_tot = ev(S['S_tot']); v_loc = ev(S['S_local']); v_s1=ev(S['S_s1']); v_s2=ev(S['S_s2'])
print("eval %.1fs"%(__import__('time').time()-t0))
P3 = mp.mpf(pickle.load(open('oct5/P3_720.pkl','rb'))['P3'])
mixed_target = mp.mpf(pickle.load(open('oct5/P3_720.pkl','rb'))['mixed'])
exLS_target  = mp.mpf(pickle.load(open('oct5/P3_720.pkl','rb'))['exact_LS'])
print("mixed:  mine  =", mp.nstr(um*v_mix, 40))
print("        target=", mp.nstr(mixed_target, 40))
print("        |diff| =", mp.nstr(abs(um*v_mix-mixed_target),3))
print("loc+s1+s2 (u*):", mp.nstr(um*(v_loc+v_s1+v_s2), 40))
print("exact_LS target:", mp.nstr(exLS_target,40))
print("        |diff| =", mp.nstr(abs(um*(v_loc+v_s1+v_s2)-exLS_target),3))
print("P3 mine = 48*u*S_tot =", mp.nstr(48*um*v_tot, 40))
print("P3 720             =", mp.nstr(P3,40))
print("|diff| =", mp.nstr(abs(48*um*v_tot-P3),3))
with open('progress.log','a') as f:
    f.write("step_03/04: mixed_engine.py exact IBP engine (20 mixed integrals verified 103 digits vs quadrature);\n"
            "  assembled S_tot (2622 terms); mixed u*S_mix vs target diff=%s; P3=48u*S_tot vs P3_720 diff=%s\n"
            % (mp.nstr(abs(um*v_mix-mixed_target),3), mp.nstr(abs(48*um*v_tot-P3),3)))
