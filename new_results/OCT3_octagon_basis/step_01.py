import os, json
from mpmath import mp, mpf, sqrt, sin, cos, pi, nstr

mp.dps = 60
print("CWD:", os.getcwd())
print("Files here:", os.listdir('.'))

path = 'research/08_theory/D15b_polygon_law/vals_m2_8_200.json'
print("Data file exists:", os.path.exists(path))

if os.path.exists(path):
    data = json.load(open(path))
    print("type:", type(data), "keys:", list(data.keys())[:10])
    v = data["8"]
    print("entry type:", type(v), "len:", len(v))
    for i, e in enumerate(v):
        print(i, repr(e)[:100])

# Identification checks
R = sqrt(2) - 1                      # circumradius of regular octagon K_1
P1_oct = 8 * R * sqrt(2 - sqrt(2))   # perimeter = 8 * side, side = R*sqrt(2-sqrt(2))
print("P(1) [regular octagon, circumradius sqrt2-1] =", nstr(P1_oct, 50))

C4 = 192 * sin(pi/8)**4 * cos(pi/8) * (3*sin(pi/8)**2 - 1)
print("C4 (local edge law zeta(3) coefficient) =", nstr(C4, 60))
print("C4 simplified = 6*b^2*cos(pi/8)*(2-3*sqrt2):")
b = sqrt(2)-1
print(nstr(6*b**2*cos(pi/8)*(2-3*sqrt(2)), 60))

target = "-4.39449630527605215997911011380609991382237744601256533876135945109576"
print("target P'''(1) prefix:", target)
