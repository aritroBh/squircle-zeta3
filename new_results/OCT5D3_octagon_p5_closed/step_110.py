import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
from fractions import Fraction as F
import mpmath as mp
from oct5d3_atoms import prod

FORMS = {}
FORMS['Wt1'] = {'lam': ('4800', '-3360'), 'L': ('2880', '-1920')}
FORMS['Wt2'] = {'pi2': ('360', '-280'), 'lam^2': ('-20880', '14760'), 'lam*L': ('-19200', '13440'),
                'L^2': ('-7680', '5280'), 'Li2(b)': ('-3840', '2880')}
FORMS['Wt3'] = {'Li3(b)': ('-4320', '2880'), 'Li3(-b)': ('15120', '-10800'), 'Li3(1/sqrt2)': ('11520', '-8640'),
                'z3': ('-17985/2', '6915'), 'lam^0*L^3': ('6480', '-4560'), 'lam^1*L^2': ('17280', '-12240'),
                'lam^2*L^1': ('34200', '-24120'), 'lam^3*L^0': ('34410', '-24330'),
                'Li2(b)*lam': ('5760', '-4320'), 'Li2(b)*L': ('5760', '-4320'),
                'pi2*lam': ('1370', '-980'), 'pi2*L': ('-900', '660')}
FORMS['Wt4'] = {'Li4(b)': ('56880', '-38520'), 'Li4(-b)': ('41880', '-27720'), 'Li4(-b^2)': ('-5100', '3450'),
                'Li4(1/sqrt2)': ('40920', '-29400'), 'Li4(1/2)': ('27930', '-19680'), 'Li4(1-1/sqrt2)': ('9240', '-6600'),
                'Li4(2-sqrt2)': ('-7200', '4800'), 'lam^0*L^4': ('-1820', '1255'), 'lam^1*L^3': ('-3520', '2620'),
                'lam^2*L^2': ('-11010', '7755'), 'lam^3*L^1': ('-66175/2', '46745/2'), 'lam^4*L^0': ('-67825/2', '191835/8'),
                'pi2*lam^0*L^2': ('1080', '-1555/2'), 'pi2*lam^1*L^1': ('-2020', '1405'), 'pi2*lam^2*L^0': ('-5305/4', '950'),
                'pi2^2': ('-18421/24', '3311/6'), 'Li2(b)^2': ('-1080', '600'), 'Li2(b)*lam^0*L^2': ('-3000', '2040'),
                'Li2(b)*lam^1*L^1': ('-4320', '3240'), 'Li2(b)*lam^2*L^0': ('3360', '-2310'), 'pi2*Li2(b)': ('2060', '-1465'),
                'z3*lam': ('64665/4', '-45705/4'), 'z3*L': ('7395', '-5610'), 'Li3(b)*lam': ('-7440', '5040'),
                'Li3(b)*L': ('2880', '-1920'), 'Li3(-b)*lam': ('-3540', '2580'), 'Li3(-b)*L': ('-11760', '8400'),
                'Li3(1/sqrt2)*lam': ('15840', '-11040'), 'Li3(1/sqrt2)*L': ('-8640', '6480')}
pickle.dump(FORMS, open("closed_forms_lower.pkl","wb"))

# specs for evaluation: map atom label -> spec (reuse W-basis specs structure)
def spec_of(label):
    # labels like 'lam^2*L^1', 'pi2*lam^2*L^0', 'Li3(b)*lam', 'z3*L', 'pi2^2', 'Li2(b)^2'
    parts = label.split('*')
    out = []
    for p in parts:
        if p.startswith('lam^'): out.append(('lam', int(p[4:])))
        elif p == 'lam': out.append(('lam', 1))
        elif p.startswith('L^'): out.append(('L', int(p[2:])))
        elif p == 'L': out.append(('L', 1))
        elif p == 'pi2': out.append(('pi2', 1))
        elif p == 'pi2^2': out.append(('pi2', 2))
        elif p == 'z3': out.append(('z3', 1))
        elif p.startswith('Li2(b)^'): out.append(('Li2(b)', int(p[7:])))
        else: out.append((p, 1))
    return out

for D in (1300,):
    mp.mp.dps = D
    A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
    Wt = pickle.load(open("wparts_true_1650.pkl","rb"))
    r2 = mp.sqrt(2)
    for lvl in (1, 2, 3, 4):
        tot = mp.mpf(0)
        for label, (a, b) in FORMS[f'Wt{lvl}'].items():
            tot += (mp.mpf(F(a).numerator)/F(a).denominator + mp.mpf(F(b).numerator)/F(b).denominator*r2)*prod(A, spec_of(label))
        res = abs(tot - mp.mpf(Wt[f'Wt{lvl}']))
        print(f"Wt{lvl} direct check at {D}: residual", mp.nstr(res, 3))
print("saved closed_forms_lower.pkl")
