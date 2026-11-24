import sympy as s
from itertools import product, combinations_with_replacement, permutations

t = s.symbols('t', positive=True)
letters = [t, 1-t, 1+t, s.Integer(2)]

def logs(f):
    out = [0]*4
    num, den = s.fraction(s.factor(f))
    for a, sgn in [(num, 1), (den, -1)]:
        const, facs = s.factor_list(a, t)
        const = abs(const)
        out[3] += 0*sgn*(s.factorint(s.numer(const)).get(2, 0)-s.factorint(s.denom(const)).get(2, 0))
        for fac, p in facs:
            j = next((j for j, z in enumerate(letters[:3]) if s.simplify(fac/z).is_number), None)
            if j is None:
                raise ValueError((f, fac))
            out[j] += sgn*p
    return s.Matrix(out)

def tens(*v):
    return s.Matrix([s.prod(v[j][i] for j, i in enumerate(q)) for q in product(range(4), repeat=len(v))])

def sy3(z):
    return -tens(logs(1-z), logs(z), logs(z))

def prod21(z, l):
    a, b, c = logs(1-z), logs(z), logs(l)
    return -tens(a,b,c)-tens(a,c,b)-tens(c,a,b)

def prod111(*l):
    v = list(map(logs,l))
    return sum((tens(*(v[j] for j in p)) for p in permutations(range(3))), s.zeros(64,1))

R,A,B,C = map(logs, letters)
c = (1-t)/(1+t)
b = (1-t)/2
J = tens(B-C-R,B)+tens(B,B-C-R)-tens(C+R-B,A-B)+tens(R,A)
JH = -tens(A-B,B-C)-tens(B-C,A-B)-tens(A-B,C+R-B)+tens(A,R)
F = s.Matrix([JH[i]*A[j] for i in range(16) for j in range(4)])-tens(B-C,A-C,A-B)+tens(B,R,R-B)
zlist = [t,-t,c,-c,1-t,b,(1+t)/2,1/(1+t),t/(1+t),2*t/(1+t),-(1-t)/(2*t),-2*t/(1-t)]
products = [(z,l) for z in zlist for l in letters]
pure = list(combinations_with_replacement(letters,3))
cols = [sy3(z) for z in zlist]+[prod21(z,l) for z,l in products]+[prod111(*l) for l in pure]
labels = [s.polylog(3,z) for z in zlist]+[s.polylog(2,z)*s.log(l) for z,l in products]+[s.prod(s.log(v) for v in l) for l in pure]
M = s.Matrix.hstack(*cols)
print('symbol rank =', M.rank(), 'augmented rank =', M.row_join(F).rank())
sol, params = M.gauss_jordan_solve(F)
sol = sol.subs({p:0 for p in params})
E = 0
for n, v in zip(labels, sol):
    if v:
        print(v, n)
        E += v*n
print('E =', E)

# Differentiate using real logarithms, then eliminate dilogarithms by
# differentiated Euler, Landen and the c=(1-t)/(1+t) identity.
R0,A0,B0,C0,Z = s.symbols('R A B C Z')
Dt,Dm,Db = s.symbols('Dt Dm Db')
def logexact(f):
    num,den = s.fraction(s.factor(f))
    value = 0
    for part,sgn in ((num,1),(den,-1)):
        const,facs = s.factor_list(part,t)
        const = abs(const)
        value += sgn*s.log(const)/s.log(2)*C0
        for fac,p in facs:
            j = next(j for j,z in enumerate(letters[:3]) if s.simplify(fac/z).is_number)
            value += sgn*p*[R0,A0,B0][j]
    return value
Dc = -Dt+Dm-Db-R0*A0+R0*B0-B0**2/2+C0*B0-C0**2/2+3*Z/2
Dmap = {
    t:Dt, -t:Dm, c:Dc, -c:-Db-(C0-B0)**2/2,
    1-t:Z-Dt-R0*A0,
    b:Db, (1+t)/2:Z-Db-(B0-C0)*(A0-C0),
    1/(1+t):Z+Dm+R0*B0-B0**2/2,
    t/(1+t):-Dm-B0**2/2,
    2*t/(1+t):Z-Dc-(C0+R0-B0)*(A0-B0),
    -(1-t)/(2*t):-Dc-(C0+R0-B0)**2/2,
    -2*t/(1-t):-Dmap_placeholder if False else -(Z-Dc-(C0+R0-B0)*(A0-B0))-(B0-A0)**2/2,
}
Tvars = s.symbols('T:'+str(len(zlist)))
Dvars = s.symbols('D:'+str(len(zlist)))
formal = E.xreplace({**{s.polylog(3,z):v for z,v in zip(zlist,Tvars)},
                     **{s.polylog(2,z):v for z,v in zip(zlist,Dvars)},
                     **dict(zip(map(s.log,letters),[R0,A0,B0,C0]))})
DE = sum(s.diff(formal,Tv)*Dv*s.diff(z,t)/z-s.diff(formal,Dv)*logexact(1-z)*s.diff(z,t)/z for z,Tv,Dv in zip(zlist,Tvars,Dvars))
DE += sum(s.diff(formal,v)*s.diff(l,t)/l for l,v in zip(letters,[R0,A0,B0,C0]))
DE = DE.subs(dict(zip(Dvars,[Dmap[z] for z in zlist])))
JG = -(A0-B0)*(B0-C0)+Dmap[2*t/(1+t)]-Dt
DG = -JG/(1-t)-2*Db/(1-t*t)-Dm/(t*(1+t))
res = s.factor(DE-DG)
print('derivative E-H =',res)
H0 = C0*B0**2/2-C0**2*B0/2+Z*A0/4+3*Z*B0/4
dH0 = sum(s.diff(H0,v)*s.diff(l,t)/l for l,v in zip(letters,[R0,A0,B0,C0]))
print('residual after H0 =',s.factor(res-dH0))
P = formal.subs(dict(zip(Tvars,[0]*len(zlist)))).subs(dict(zip(Dvars,[Dmap[z] for z in zlist])))
Q = s.factor(P-H0)
print('H product Q =',s.collect(Q,[Dt,Dm,Db,Z]))
Bk,Rk = s.symbols('Bk Rk')
print('Q at octagon logs =',s.collect(s.expand(Q.subs({A0:Bk+Rk,B0:Bk,C0:2*Bk,R0:Rk})),[Dt,Dm,Db,Z]))

if __name__ == '__main__':
    import mpmath as mp
    from verify import depth2
    mp.mp.dps = 80
    Efun = s.lambdify(t, E, 'mpmath')
    for u in map(mp.mpf, ['0.2','0.4','0.6','0.8']):
        f = depth2((1-u)/2,-2*u/(1-u))
        print('G-E at',u,'=',mp.nstr(f-Efun(u),40))
