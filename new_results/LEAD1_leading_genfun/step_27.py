import mpmath, json
mpmath.mp.dps = 60
from mpmath import mpf, sqrt, log, quad, pi, sin, tan, cos
a = 1+sqrt(2); b = sqrt(2)-1
def t_of_v(v): return (v-b)*(v+a)/(4*v)
def u_of_v(v): return (a-v)*(v+b)/(4*v)
def v_int(x,y,R):
    return quad(lambda v: t_of_v(v)**x * u_of_v(v)**y * R(v), [b,1,a])
def psi_int(x,y,M):
    return quad(lambda ps: ((1+tan(ps))/2)**x * ((1-tan(ps))/2)**y * M(ps), [-pi/4, 0, pi/4])
combos = {
 'odd: N_-a - N_a': (lambda v: 1/(v+a)-1/(v-a), lambda ps: 1/(1-sqrt(2)*sin(ps))),
 'even: N_-a - N_-b': (lambda v: 1/(v+a)-1/(v+b), None),
 'diag: N_-a - N_0/2': (lambda v: 1/(v+a)-1/(2*v), None),
 'oddlast: N_-a - N_a/2 + N_-b/2 - N_0/2': (lambda v: 1/(v+a)-1/(2*(v-a))+1/(2*(v+b))-1/(2*v), None),
}
# first: verify odd combo with corrected factor
for x,y in [(1,0),(2,1),(0.5,0.25),(3,0)]:
    d = v_int(x,y,combos['odd: N_-a - N_a'][0]) - psi_int(x,y,combos['odd: N_-a - N_a'][1])
    print('odd', x,y, mpmath.nstr(d,8))
