import mpmath, json
mpmath.mp.dps = 60
from mpmath import mpf, sqrt, log, quad, pi, sin, tan, cos
os_import = __import__('os'); os_import.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
a = 1+sqrt(2); b = sqrt(2)-1
def t_of_v(v): return (v-b)*(v+a)/(4*v)
def u_of_v(v): return (a-v)*(v+b)/(4*v)
# measure identity test: odd d combo N_{-a}-N_{a} should equal int t^x(1-t)^y * 2/(1-sqrt2 sin psi) dpsi
def lhs(x,y):  # int_b^a t^x(1-t)^y [1/(v+a)-1/(v-a)] dv
    return quad(lambda v: t_of_v(v)**x * u_of_v(v)**y * (1/(v+a)-1/(v-a)), [b,1,a])
def rhs(x,y):
    return quad(lambda ps: ((1+tan(ps))/2)**x * ((1-tan(ps))/2)**y * 2/(1-sqrt(2)*sin(ps)), [-pi/4, 0, pi/4])
for x,y in [(1,0),(2,0),(1,1),(3,0),(2,1),(0.5,0.25)]:
    print(x,y, mpmath.nstr(lhs(x,y)-rhs(x,y), 8))
