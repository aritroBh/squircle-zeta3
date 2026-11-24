import mpmath as mp
mp.mp.dps = 50
x = mp.mpf('0.2'); y = mp.mpf('0.3'); xi = 1-x; eta = 1-y
L3 = lambda z: mp.polylog(3, z); L4 = lambda z: mp.polylog(4, z)
K20 = [
    (-1,2,1,-1,1, 1), (-1,1,2,1,-1, 1), (1,2,1,-1,-2, 1), (1,1,2,2,-1, 1),
    (1,1,1,0,0, -6), (1,1,1,-1,-1, -6), (-1,1,1,-1,0, -6), (-1,1,1,0,-1, -6),
    (1,1,0,0,1, -3), (1,0,1,1,0, -3), (-1,1,0,-1,1, -3), (-1,0,1,1,-1, -3),
    (-1,1,0,-1,-1, -3), (-1,0,1,-1,-1, -3), (1,1,0,0,-1, -3), (1,0,1,-1,0, -3),
    (1,1,0,0,0, 6), (1,0,1,0,0, 6), (-1,1,0,-1,0, 6), (-1,0,1,0,-1, 6),
]
def mon(s,a,b,c,d): return s*x**a*y**b*xi**c*eta**d
# K4 itself
G = sum(di*L4(mon(s,a,b,c,d)) for s,a,b,c,d,di in K20) - mp.mpf('1.5')*mp.log(xi)**2*mp.log(eta)**2
print('K4 residual:', mp.nstr(G, 6))
# derivative wrt x
Rx = sum(di*(a*xi - c*x)*L3(mon(s,a,b,c,d)) for s,a,b,c,d,di in K20) + 3*x*mp.log(xi)*mp.log(eta)**2
print('Rx residual:', mp.nstr(Rx, 6))
# derivative wrt y
Ry = sum(di*(b*eta - d*y)*L3(mon(s,a,b,c,d)) for s,a,b,c,d,di in K20) + 3*y*mp.log(xi)**2*mp.log(eta)
print('Ry residual:', mp.nstr(Ry, 6))
# check derivative numerically
h = mp.mpf('1e-20')
Gx = (sum(di*L4((lambda xx,xix: s*xx**a*y**b*xix**c*eta**d)(x+h, 1-(x+h))) for s,a,b,c,d,di in K20)
      - mp.mpf('1.5')*mp.log(1-(x+h))**2*mp.log(eta)**2 - 
      (sum(di*L4((lambda xx,xix: s*xx**a*y**b*xix**c*eta**d)(x-h, 1-(x-h))) for s,a,b,c,d,di in K20)
      - mp.mpf('1.5')*mp.log(1-(x-h))**2*mp.log(eta)**2))/(2*h)
print('numeric dG/dx:', mp.nstr(Gx, 6))
