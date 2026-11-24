import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import mpmath as mp
mp.mp.dps = 250
r2 = mp.sqrt(2); X = mp.log(2); Y = mp.log(1+r2); b = r2-1
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
av = pickle.load(open("atom_vals.pkl","rb"))
Li5v = {k: mp.mpf(v) for k,v in av['Li5v'].items()}
Li4v = {k: mp.mpf(v) for k,v in av['Li4v'].items()}
Li3v = {k: mp.mpf(v) for k,v in av['Li3v'].items()}
Li2v = {k: mp.mpf(v) for k,v in av['Li2v'].items()}
L32_b1 = mp.mpf(av['L32_b1']); L32_h1 = mp.mpf(av['L32_h1'])

data = pickle.load(open("symbol_fit2.pkl","rb"))
cols = data['cols']
xsol = pickle.load(open("symbol_fit2_solution.pkl","rb"))['x']
nz = [(label, kind, xv) for (label, kind), xv in zip(cols, xsol) if xv.p != 0 or xv.q != 0]
import re
def term_value(label):
    m = re.match(r"Li5\((.*)\)$", label)
    if m: return Li5v[m.group(1)]
    m = re.match(r"Li4\((.*)\)\*([XY])$", label)
    if m: return Li4v[m.group(1)] * (X if m.group(2)=='X' else Y)
    for pat, fac in [(r"Li3\((.*)\)\*X\^2/2$", X**2/2), (r"Li3\((.*)\)\*X\*Y$", X*Y), (r"Li3\((.*)\)\*Y\^2/2$", Y**2/2)]:
        m = re.match(pat, label)
        if m: return Li3v[m.group(1)]*fac
    for pat, fac in [(r"Li2\((.*)\)\*X\^3/6$", X**3/6), (r"Li2\((.*)\)\*X\^2Y/2$", X**2*Y/2),
                     (r"Li2\((.*)\)\*XY\^2/2$", X*Y**2/2), (r"Li2\((.*)\)\*Y\^3/6$", Y**3/6)]:
        m = re.match(pat, label)
        if m: return Li2v[m.group(1)]*fac
    m = re.match(r"Li3\((.*)\)\*Li2\((.*)\)$", label)
    if m: return Li3v[m.group(1)] * Li2v[m.group(2)]
    m = re.match(r"Li2\((.*)\)\*Li2\((.*)\)\*([XY])$", label)
    if m: return Li2v[m.group(1)] * Li2v[m.group(2)] * (X if m.group(3)=='X' else Y)
    if label == "Li32(b,1)": return L32_b1
    if label == "Li32(1/sqrt2,1)": return L32_h1
    raise ValueError(label)

cand = mp.mpf(0)
print("%-32s %-22s %s" % ("label","coeff","coeff*value"))
for label, kind, xv in nz:
    tv = term_value(label)
    cv = q2f(xv)*tv
    cand += cv
    print("%-32s %-22s %s" % (label, str(xv), mp.nstr(cv, 12)))
W5 = mp.mpf(pickle.load(open("w5_250.pkl","rb"))['W5'])
print("candidate =", mp.nstr(cand, 20))
print("W5        =", mp.nstr(W5, 20))
print("R         =", mp.nstr(W5-cand, 20))
