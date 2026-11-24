import sys, os, time, pickle
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 4200
cc.N = 1400
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 60

I = acb(0,1)
two_pi_i = 2*arb.pi()*I

def letter(eta_sign):  # eta = sign*xi8 -> engine letter 1 (+xi8) or 5 (-xi8)
    return 1 if eta_sign==1 else 5

def zeta_word(xs, eps):
    """Glanois zeta(x1..xp ; eps1..eps_{p-1}, eps_p xi8): word for engine h().
    Convention (1.1): sum over 0<n1<...<np of prod eps_k^{n_k}/n_k^{x_k}."""
    p = len(xs)
    z = list(eps[:-1]) + [eps[-1]]  # eps_k = +-1 for k<p, last factor eps_p multiplies xi8
    # eta_k = (eps_k*...*eps_p)*xi8 ; sign = prod_{j>=k} eps_j
    T = []
    prod = 1
    signs = [1]*(p)
    run = 1
    for k in range(p-1,-1,-1):
        run *= z[k]
        signs[k] = run
    for k in range(p):
        T.append(letter(signs[k]))
        T.extend([-1]*(xs[k]-1))
    return tuple(reversed(T))

# Build the 27 basis elements of B^8 at weight 3:  wt = sum xi + s  (wt(2 pi i)=1)
import itertools
basis = []  # (label, xs, eps, s)
for p in (1,2,3):
    for xs in itertools.product([1,2,3], repeat=p):
        if sum(xs) != 3: continue
        for eps in itertools.product([1,-1], repeat=p):
            basis.append((f'z{xs};{eps}', xs, eps, 0))
for p in (1,2):
    for xs in itertools.product([1,2], repeat=p):
        if sum(xs) != 2: continue
        for eps in itertools.product([1,-1], repeat=p):
            basis.append((f'z{xs};{eps}*T', xs, eps, 1))
for eps in itertools.product([1,-1], repeat=1):
    basis.append((f'z(1);{eps}*T^2', (1,), eps, 2))
basis.append(('T^3', (), (), 3))
print('basis size:', len(basis))
assert len(basis)==27

t=time.time()
basis_vals = []
basis_words = []
for (label, xs, eps, s) in basis:
    if xs:
        w = zeta_word(xs, eps)
        v = cc.h(w)
    else:
        w = ()
        v = acb(1)
    v = v * two_pi_i**s
    basis_words.append(w); basis_vals.append(v)
    print(label, w, 'rad', v.real.rad().str(2,radius=False))
print('eval time', time.time()-t)

# store at 600 digits
def s600(z):
    return (z.real.mid().str(600,radius=False,more=True), z.imag.mid().str(600,radius=False,more=True))
with open(WORK+'/basis_vals.pickle','wb') as f:
    pickle.dump({'basis':[b[0] for b in basis],'words':basis_words,
                 'vals':[s600(v) for v in basis_vals],
                 'info':{'ctx.prec':4200,'N':1400}},f)
# quick external validation: shuffle product check  z(1;xi)z(2;xi) = z(1,2;1,xi)+z(2,1;1,xi)+z(3;xi)
ix = {b[0]:i for i,b in enumerate(basis)}
L = basis_vals[ix["z(1,);(1,)"]]*basis_vals[ix["z(2,);(1,)"]]
R = basis_vals[ix["z((1, 2));(1, 1)"]]+basis_vals[ix["z((2, 1));(1, 1)"]]+basis_vals[ix["z(3,);(1,)"]]
d = L-R
print('stuffle1 residual:', mp.nstr(abs(cc.approx(d)),3))
with open(WORK+'/progress.log','a') as f:
    f.write('step_05: built 27-element B^8 basis at weight 3, evaluated at prec 4200 bits/N=1400, saved basis_vals.pickle; stuffle check\n')
