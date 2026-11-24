import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
from fractions import Fraction as F
from oct5d3_basis import atom_specs, C_SPECS
from oct5d3_atoms import INVISIBLE5, prod, k32_value

sf = pickle.load(open("symfit_m4f.pkl","rb"))
x0, labels = sf['x0'], sf['vis_labels']
Kred = pickle.load(open("kernel_red.pkl","rb"))['Kred']
mixed = pickle.load(open("mixed61_vals.pkl","rb"))
spec_map = {lab: spec for lab, kind, spec in atom_specs()}
vis = [(lab, spec_map[lab]) for lab in labels]
def q2f(z):
    return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*mp.sqrt(2)

# the 1500-digit relation (from previous cell), hard-coded
REL = {
 'z2*Li3(b)': (F(-1616917,184), F(172877,46)), 'z2*Li3(-b)': (F(-6566517,184), F(568116,23)),
 'z2*Li3(1/sqrt2)': (F(-415,4), F(1735,2)), 'z2*z3': (F(-16159155,2944), F(-8342511,2944)),
 'z2*lam^3': (F(-359551,12), F(2057887,96)), 'z2*lam^2*L': (F(2321073,32), F(-807141,16)),
 'z2*lam*L^2': (F(-3581371,32), F(78827,1)), 'z2*L^3': (F(39019803,736), F(-10222387,276)),
 'z2*Li2(b)*lam': (F(-28305,2), F(40005,4)), 'z2*Li2(b)*L': (F(-2520,1), F(1845,1)),
 'z2*pi2*lam': (F(91133485,8832), F(-30695753,4416)), 'z2*pi2*L': (F(-26007269,8832), F(812637,368)),
 'z3*pi2': (F(7470059,2944), F(6472285,8832)), 'z3*lam^2': (F(189115,128), F(-26575,16)),
 'z3*lam*L': (F(542015,32), F(-514945,32)), 'z3*L^2': (F(-4742397,368), F(6671861,736)),
 'z3*Li2(b)': (F(-66645,8), F(93195,16)), 'z4*lam': (F(-698145,92), F(-1364255,2944)),
 'z4*L': (F(6672697,1472), F(-10009091,2944)), 'z5': (F(82937921,2944), F(-56662239,2944)),
 'Lchi5': (F(-400,3), F(39752,3)), 'K32': (F(-3272927,1104), F(7011289,2944)),
 'MIX1': (F(-4507895,138), F(2395121,92)), 'MIX2': (F(2941837,69), F(-6339427,184)),
 'MIX3': (F(2032204,69), F(-3665635,184)), 'MIX4': (F(10695,1), F(-10980,1)),
 'MIX5': (F(28095,2), F(-22545,2)), 'MIX6': (F(-16450,1), F(10610,1)),
 'MIX7': (F(-1417693,1104), F(2744391,2944)),
 'K1': (F(-95830013,384), F(64714967,368)), 'K2': (F(3431650,69), F(-92894123,2944)),
 'K3': (F(-281424457,4416), F(4514659,92)), 'K4': (F(-10096565,2944), F(13035935,8832)),
 'K5': (F(-300,1), F(210,1)), 'K6': (F(-1656531,184), F(1395579,184)),
 'K7': (F(28780351,2944), F(-15206245,2208)), 'K8': (F(-21877223,2208), F(3024617,384)),
 'K9': (F(-31884357,2944), F(7588411,8832)), 'K10': (F(104,23), F(-74,23)),
 'K11': (F(-6215,1472), F(14075,4416)), 'K12': (F(-16175,736), F(12413,736)),
 'K13': (F(-456541,8832), F(6507,184)), 'K14': (F(19697995,8832), F(-2328855,1472)),
 'K15': (F(-284063,1104), F(538037,2944)), 'K16': (F(12117,184), F(-9321,184)),
 'K17': (F(-383599,8832), F(273535,8832)), 'K18': (F(15209915,8832), F(-2691331,2208)),
 'K19': (F(-94117,552), F(357871,2944)), 'K20': (F(-1234501,8832), F(440545,4416)),
 'K21': (F(-1681,184), F(1089,184)), 'K22': (F(6257,4416), F(-359,368)),
}

def setup(D):
    mp.mp.dps = D
    A = {k: mp.mpf(v) for k, v in pickle.load(open("atoms_1650.pkl","rb")).items()}
    V = [prod(A, spec) for lab, spec in vis]
    Wt5 = mp.mpf(pickle.load(open("wparts_true_1650.pkl","rb"))['Wt5'])
    classical = sum(q2f(x0[j])*V[j] for j in range(len(V)))
    Cval = sum(q2f(cf)*prod(A, spec) for lab, spec, cf in C_SPECS)
    R = Wt5 - classical - Cval
    r2 = mp.sqrt(2)
    kvals = []
    for row in Kred:
        v = mp.mpf(0)
        for j in range(len(labels)):
            p, q = row[2*j], row[2*j+1]
            if p or q: v += (mp.mpf(p) + mp.mpf(q)*r2)*V[j]
        kvals.append(v)
    vals = dict([(nm, prod(A, sp)) for nm, sp in INVISIBLE5]
                + [("K32", k32_value())]
                + [(f"MIX{i+1}", mp.mpf(val)) for i, (val, desc) in enumerate(mixed)]
                + [(f"K{i+1}", v) for i, v in enumerate(kvals)])
    return R, vals

for D in (900, 1500):
    R, vals = setup(D)
    r2 = mp.sqrt(2)
    tot = mp.mpf(0)
    for nm, (p, q) in REL.items():
        tot += (mp.mpf(p.numerator)/p.denominator + mp.mpf(q.numerator)/q.denominator*r2)*vals[nm]
    print(f"direct check at {D}: |R - relation| =", mp.nstr(abs(R-tot), 3))
