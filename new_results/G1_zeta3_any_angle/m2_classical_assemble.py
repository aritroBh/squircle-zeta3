"""m2_classical_assemble.py -- the uniform CLASSICAL formula for M2(alpha), all alpha in (0, pi).

  M2 = z2*Phi(c) + R2(c),  R2 = sum_w p_w(c) H_w(c)   (hpl_fit_k2_(0,0,0).json, p_w in Q[c, 1/c]),
  each weight-3 H_w(x) replaced by its classical form (classical_dict.json: Li3 at the 8 arguments
  x, -x, 1-x, 1/(1+x), (1+x)/2, (1-x)/(1+x), -(1-x)/(1+x), 2x/(1+x), Li2*log, log^3, pi^2 log, zeta(3));
  weight 2: H_{0,1} = Li2(x), H_{0,-1} = -Li2(-x).
Outputs
 (1) the coefficient of each classical function as an element of Q[c,1/c] (sympy), in particular the
     coefficient Z(c) of zeta(3) in this uniform classical basis;
 (2) a 300-digit check of the classical formula against the stored M2 at the 9 angles;
 (3) the explicit closed forms at pi/3 (c = sqrt3/2) and 2pi/3 (c = 1/2), coefficients in Q(c).
"""
import json, re
import sympy as sp
import mpmath as mp

c, x = sp.symbols('c x', positive=True)
R2 = json.load(open('hpl_fit_k2_(0,0,0).json'))
CD = json.load(open('classical_dict.json'))
CONSTS = {'1': sp.Integer(1), 'log2': sp.log(2), 'pi2': sp.pi**2, 'log2^2': sp.log(2)**2, 'z3': sp.zeta(3),
          'pi2*log2': sp.pi**2*sp.log(2), 'log2^3': sp.log(2)**3}
ARGS = {'x': x, '-x': -x, '1-x': 1 - x, '1/(1+x)': 1/(1 + x), '(1+x)/2': (1 + x)/2,
        '(1-x)/(1+x)': (1 - x)/(1 + x), '-(1-x)/(1+x)': -(1 - x)/(1 + x), '2x/(1+x)': 2*x/(1 + x)}
FUN = {'lx': sp.log(x), 'l1m': sp.log(1 - x), 'l1p': sp.log(1 + x)}


def fun_expr(name):
    """sympy expression (in x) of a classical basis function name from classical_dict."""
    if name == '1':
        return sp.Integer(1)
    parts = name.split('*')
    out = sp.Integer(1)
    for p in parts:
        m = re.match(r'Li([23])\((.*)\)', p)
        if m:
            out *= sp.polylog(int(m.group(1)), ARGS[m.group(2)])
        else:
            out *= FUN[p]
    return out


def arg_of(name):
    m = re.match(r'Li([23])\((.*)\)', name)
    return m.group(2) if m else None


# --- collect: key = (function name in x), value = coefficient in Q[c,1/c] * constant
coef = {}
def add(fname, cconst, val):
    key = (fname, cconst)
    coef[key] = sp.simplify(coef.get(key, 0) + val)

for name, cstr in R2:
    mm = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    e = int(mm.group(1)); w = [int(z) for z in mm.group(2).split(',')]
    pref = sp.Rational(cstr)*c**e
    if w == [0, 1]:
        add('Li2(x)', '1', pref)
    elif w == [0, -1]:
        add('Li2(-x)', '1', -pref)
    else:
        for fn, cn, fr in CD[str(w)]:
            add(fn, cn, pref*sp.Rational(fr))

print('UNIFORM CLASSICAL FORM:  M2 = z2*Phi(c) + sum  coef(c) * const * F(x=c)')
for (fn, cn), v in sorted(coef.items(), key=lambda kv: (kv[0][0], kv[0][1])):
    if v != 0:
        print(f'   [{cn:>9}] * {fn:<22} : {sp.factor(v)}')
Z = sp.factor(coef.get(('1', 'z3'), 0))
print('\nCoefficient of zeta(3) in this uniform classical basis:  Z(c) =', Z, '=', sp.expand(Z))

# --- numerical check at 300 digits
mp.mp.dps = 320
num_terms = []
for (fn, cn), v in coef.items():
    if v == 0:
        continue
    fx = sp.lambdify(x, fun_expr(fn), 'mpmath')
    cv = sp.lambdify(c, v, 'mpmath')
    kv = {'1': 1, 'log2': mp.log(2), 'pi2': mp.pi**2, 'log2^2': mp.log(2)**2, 'z3': mp.zeta(3),
          'pi2*log2': mp.pi**2*mp.log(2), 'log2^3': mp.log(2)**3}[cn]
    num_terms.append((fx, cv, kv))


def M2_classical(cc):
    Phi = mp.mpf(1)/2 - (1 - cc*cc)/(2*cc)*mp.atanh(cc)
    return mp.zeta(2)*Phi + mp.fsum(cv(cc)*kv*fx(cc) for fx, cv, kv in num_terms)


ANG = {'pi_2': mp.pi/2, '2pi_3': 2*mp.pi/3, 'pi_3': mp.pi/3, '3pi_4': 3*mp.pi/4, 'pi_6': mp.pi/6,
       '5pi_6': 5*mp.pi/6, 'pi_4': mp.pi/4, 'gen1p3': mp.mpf(13)/10, 'gen2p0': mp.mpf(2)}
for tag, a in ANG.items():
    st = [l for l in open(f'values_{tag}.txt') if l.startswith('M2 =')][0]
    ref = mp.mpf(st.split('=')[1].strip())
    cc = mp.cos(a/2)
    zval = sp.lambdify(c, Z, 'mpmath')(cc)
    print(f'{tag:>7}: |classical formula - series M2| = {mp.nstr(abs(M2_classical(cc) - ref), 3)}'
          f'   Z(c) = {mp.nstr(zval, 20)}')

# --- explicit closed forms at pi/3 and 2pi/3
for label, cval in [('2pi/3', sp.Rational(1, 2)), ('pi/3', sp.sqrt(3)/2), ('pi/2', 1/sp.sqrt(2))]:
    print(f'\nClosed form at alpha = {label} (c = {cval}):  M2 = z2*Phi + ...')
    phi = sp.nsimplify(sp.Rational(1, 2) - (1 - cval**2)/(2*cval)*sp.atanh(cval))
    print('   Phi(c) =', sp.simplify(sp.Rational(1, 2) - (1 - cval**2)/(2*cval)*sp.log((1 + cval)/(1 - cval))/2))
    grouped = {}
    for (fn, cn), v in coef.items():
        vv = sp.radsimp(sp.simplify(v.subs(c, cval)))
        if vv == 0:
            continue
        fsub = sp.simplify(fun_expr(fn).subs(x, cval))
        grouped.setdefault(str(fsub) + ('' if cn == '1' else '  *' + cn), 0)
        grouped[str(fsub) + ('' if cn == '1' else '  *' + cn)] += vv
    for k2, v2 in grouped.items():
        v2 = sp.radsimp(sp.nsimplify(sp.simplify(v2)))
        if v2 != 0:
            print(f'   {str(v2):>28}  *  {k2}')
