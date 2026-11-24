import sys, os
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
d, h, i = sp.symbols('d h'), sp.symbols('h'), sp.I
d = sp.symbols('d')
m = (d - h)/2
i = sp.I

# N_l from the fits (as functions of (d,m))
N1 = 4*m + (-2+6*i)*d - 4*i*d**2
N2 = (-40*i)*m + 16*m**2 + (sp.Rational(80,3)+20*i)*d + (-16+96*i)*d*m + (-80-48*i)*d**2 + (-32*i)*d**2*m + (sp.Rational(208,3)+16*i)*d**3 - 16*d**4
N3 = (-688)*m + (-672*i)*m**2 + 64*m**3 + (344-280*i)*d + (2768+672*i)*d*m + (-96+864*i)*d*m**2 + (-1384+1144*i)*d**2 + (-3264-864*i)*d**2*m + (-192*i)*d**2*m**2 + (1648-1984*i)*d**3 + (1408+192*i)*d**3*m + (-704+1600*i)*d**4 + (-192)*d**4*m + (96-544*i)*d**5 + (64*i)*d**6

for name, N, l in [('N1',N1,1),('N2',N2,2),('N3',N3,3)]:
    Nh = sp.expand(N.subs(m, (d-h)/2))
    # factor over Q(i)
    Nf = sp.factor(Nh, extension=[i])
    print(f'=== {name} in (d,h):')
    print('  expanded:', Nh)
    print('  factored:', Nf)
