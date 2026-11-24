
# step_11.py : numerical verification of bivariate circle-end coefficients b_jk via finite differences
import mpmath as mp
import sympy as sp
import pickle
mp.mp.dps = 90

def W(s1, s2):
    q = 2/(1-2*s1); u = 2/(1-2*s2)
    f = lambda th: 2*(q-1)*(mp.cos(th)**u+mp.sin(th)**u)**(1/u)*(mp.cos(th)*mp.sin(th))**(q-2)*(mp.cos(th)**q+mp.sin(th)**q)**(1/q-2)
    return mp.quad(f, [0, mp.pi/6, mp.pi/4, mp.pi/3, mp.pi/2])

# finite-difference weights, 9-point stencil, for derivatives 0..5 (8th order accurate)
import numpy as np
def fd_weights(n, h):
    # central weights for n-th derivative on stencil -4..4
    import sympy as s
    xs = list(range(-4,5))
    M = []
    for m in range(9):
        M.append([s.Rational(x)**m for x in xs])
    M = s.Matrix(M)
    rhs = [0]*9; rhs[n] = s.factorial(n)
    w = M.solve(s.Matrix(rhs))
    return [s.Rational(v) for v in w]

H = mp.mpf('0.004')
grid = {}
for i in range(-4,5):
    for j in range(-4,5):
        grid[(i,j)] = W(i*H, j*H)
weights = {n: fd_weights(n, H) for n in range(6)}
bex = {eval(k): sp.sympify(v) for k,v in pickle.load(open('scratchpad/circle_bivariate.pkl','rb')).items()}
for (j,k) in [(0,1),(0,2),(1,1),(0,3),(1,2),(2,2),(0,4),(1,3),(0,5),(2,3),(4,1)]:
    num = 0
    for i in range(-4,5):
        for l in range(-4,5):
            num += weights[j][i+4]*weights[k][l+4]*grid[(i,l)]
    num = num/H**(j+k)/mp.factorial(j)/mp.factorial(k)
    ex = mp.mpf(sp.N(bex[(j,k)], 60))
    print(f"b_({j},{k}): fd={mp.nstr(num,22)}  exact={mp.nstr(ex,22)}  diff={mp.nstr(abs(num-ex),3)}")
