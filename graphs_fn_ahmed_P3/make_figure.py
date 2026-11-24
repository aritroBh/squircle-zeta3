"""make_figure.py -- fig_fn_graph.pdf: (a) the slope-2/e point; (b) convergence of the constant C(f)."""
import json
import mpmath as mp
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fnlib import TESTS, make_callables

mp.mp.dps = 60
xs = sp.Symbol('x', positive=True)
K2 = 1 + mp.pi**2/6; K3 = 2 - mp.mpf(3)/2*mp.zeta(3); KAP = 1 - mp.log(2)
def Q2(a, b, Xi): return -((a**2 - b)/(2*a**3)*Xi**2 + b/a**3*Xi + ((a**2 - b)*K2 - a)/(2*a**3))
def Q3(a, b, c, Xi):
    return ((a**4 - 3*a**2*b - a*c + 3*b**2)/(6*a**5)*Xi**3 + (3*a**2*b + a*c - 4*b**2)/(2*a**5)*Xi**2
            + (a**4*K2 - 3*a**2*b*K2 - a*c*K2 + 3*b**2*K2 + 2*b**2)/(2*a**5)*Xi
            + (a**4*K3 + 9*a**2*b*K2 - 3*a**2*b*K3 - 6*a*b + 3*a*c*K2 - a*c*K3 - 12*b**2*K2 + 3*b**2*K3)/(6*a**5))

plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'mathtext.fontset': 'cm'})
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))

# (a)
expr = TESTS['expm1']; f, fp = make_callables(expr); n = 25
a = float(sp.N(sp.diff(expr, xs).subs(xs, 1)))
X = [i/2000 for i in range(2001)]
Y = [float(f(mp.mpf(x))**n) for x in X]
F = lambda t: (n - 1)*mp.log(f(1 - t)) + mp.log(fp(1 - t)) + mp.log(n) - mp.log(2/mp.e)
ts = float(mp.findroot(F, (mp.log(n) + KAP + mp.log(a))/(a*(n - 1))))
xstar = 1 - ts; ystar = float(f(mp.mpf(xstar))**n)
ax[0].plot([0, 1, 1], [0, 0, 1], color='0.6', ls='--', lw=0.8)
ax[0].plot(X, Y, 'k', lw=1.1, label=r'$y=f(x)^{25}$')
tt = [xstar - 0.12, xstar + 0.12]
ax[0].plot(tt, [ystar + (2/float(mp.e))*(t - xstar) for t in tt], color='C3', lw=0.9, label=r'tangent of slope $2/e$')
ax[0].plot([xstar], [ystar], 'o', color='C3', ms=3.5)
ax[0].annotate('', xy=(1, 0.55), xytext=(xstar, 0.55), arrowprops=dict(arrowstyle='<->', lw=0.7))
ax[0].text((xstar + 1)/2, 0.58, r'$1-x_n^*$', ha='center', fontsize=8)
ax[0].set_xlim(0.6, 1.02); ax[0].set_ylim(-0.02, 1.02)
ax[0].set_xlabel(r'$x$'); ax[0].set_ylabel(r'$y$')
ax[0].set_title(r'(a) $f=(e^x-1)/(e-1)$: $2-L_n\approx 1-x_n^*$', fontsize=9)
ax[0].legend(loc='upper left', fontsize=7, frameon=False)

# (b)
data = json.load(open('data_defects.json'))
names = ['x', 'x^2', 'expm1', 'x+x^3', 'mobius', 'sin4']
labels = {'x': r'$x$', 'x^2': r'$x^2$', 'expm1': r'$\frac{e^x-1}{e-1}$', 'x+x^3': r'$\frac{x+x^3}{2}$',
          'mobius': r'$\frac{2x}{1+x}$', 'sin4': r'$\frac{\sin(\pi x/4)}{\sin(\pi/4)}$'}
for k, name in enumerate(names):
    expr = TESTS[name]
    aa, bb, cc = [mp.mpf(sp.N(sp.diff(expr, xs, j).subs(xs, 1), 70)) for j in (1, 2, 3)]
    C = mp.log(aa) + KAP
    pts = sorted((int(r['n']), mp.mpf(r['D'])) for r in data if r['name'] == name)
    N1, E1, E3 = [], [], []
    for nn, D in pts:
        mu = 1/mp.mpf(nn - 1); L = mp.log(nn); Xi = L + C
        e1 = abs(aa*D/mu - L - C)
        e3 = abs((aa/mu)*(D - mu**2*Q2(aa, bb, Xi) - mu**3*Q3(aa, bb, cc, Xi)) - L - C)
        N1.append(float(nn)); E1.append(float(e1)); E3.append(float(max(e3, mp.mpf(10)**-80)))
    ax[1].loglog(N1, E1, color='C%d' % k, lw=0.9, label=labels[name])
    ax[1].loglog(N1, E3, color='C%d' % k, lw=0.9, ls='--')
ax[1].set_xlabel(r'$n$'); ax[1].set_ylabel(r'$|\widehat C_n - C(f)|$')
ax[1].set_title(r'(b) constant $C(f)=\ln f^\prime(1)+1-\ln 2$', fontsize=9)
ax[1].text(1.5e3, 1e-34, 'solid: $f\'(1)(n-1)D_n-\\ln n$\ndashed: minus $\\mu^2Q_2+\\mu^3Q_3$', fontsize=7, va='top')
ax[1].legend(fontsize=7, frameon=False, loc='lower left', ncol=3, handlelength=1.2, columnspacing=0.8)
fig.tight_layout()
fig.savefig('fig_fn_graph.pdf'); fig.savefig('fig_fn_graph.png', dpi=200)
print('xstar', xstar, 'tau*', ts)
