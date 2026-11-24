import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import sympy as sp
u = sp.symbols('u')
S2 = sp.sqrt(2)
S = sp.sqrt(1+u**2)
a = 1 + S2; b = S2 - 1
v = u + S
dvdu = sp.diff(v, u)
m_odd  = (1+S2*u/S)/(1-u**2)
m_e    = (1-S2/S)/(1-u**2)
m_diag = (1-u)/2*m_e
m_ol   = m_odd - 1/(2*(1-u))
combos = {'odd': 1/(v+a)-1/(v-a), 'even': 1/(v+a)-1/(v+b),
          'diag': 1/(v+a)-1/(2*v), 'oddlast': 1/(v+a)-1/(2*(v-a))+1/(2*(v+b))-1/(2*v)}
for name, Rv in combos.items():
    meas = sp.simplify(Rv*dvdu)
    tgt = {'odd':m_odd,'even':m_e,'diag':m_diag,'oddlast':m_ol}[name]
    diff = sp.radsimp(sp.simplify(sp.expand((meas-tgt).subs(S2, sp.Symbol('q')))))
    # numeric spot check at algebraic-independent points as backup
    import random
    ok = True
    for _ in range(5):
        uu = random.uniform(-0.9, 0.9)
        if abs(complex(sp.N(meas.subs(u,uu),30)) - abs(complex(sp.N(tgt.subs(u,uu),30))) < 1e-25: pass
        else: ok = abs(complex(sp.N((meas-tgt).subs(u,uu),30))) < 1e-25
        ok &= abs(complex(sp.N((meas-tgt).subs(u,uu),30))) < 1e-25
    print(f"{name}: symbolic-zero={diff==0}  numeric={ok}")
chk = sp.radsimp(sp.simplify(m_odd - (1/(1-u)-u*m_e)))
print("m_odd = 1/(1-u) - u m_e :", chk==0)
t_v = (v-b)*(v+a)/(4*v)
print("t(v) = (1+u)/2 :", sp.radsimp(sp.simplify(t_v-(1+u)/2))==0)
u1 = (a-v)*(v+b)/(4*v)
print("1-t = (1-u)/2 :", sp.radsimp(sp.simplify(u1-(1-u)/2))==0)
# also the other branch v_- = u - S maps to the same u; and v * v_other = -1 (deck transformation v -> -1/v)
v2 = u - S
print("v+ * v- = -1 :", sp.simplify(v*v2)==-1)
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: exact u-measures + t=(1+u)/2 + deck transformation, sympy radsimp\n")
