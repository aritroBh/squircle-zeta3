
import mpmath as mp, sys
mp.mp.dps = 80
sys.path.insert(0, "research/08_theory/D17_zeta5_polygons")
sys.path.insert(0, "research/08_theory/D15b_polygon_law")
import poly_coeffs5 as p5
import poly_coeffs as p3
print("dps:", mp.mp.dps)
for m in [2,3,4,5]:
    A = mp.pi/m
    c5 = p5.make_coeffs5(p5.regular_zeros(m), [1]*m, [0]*m)
    c3 = p3.make_coeffs(p3.regular_zeros(m), [1]*m, [0]*m)
    th = mp.mpf('0.37')*A
    v5 = c5(th); v3 = c3(th)
    res = max(abs(v5[k]-v3[k]) for k in range(4))
    print("m=%d th=%s : max |D5rec-D3old| for k=0..3 = %s"%(m, mp.nstr(th,6), mp.nstr(res,3)))
