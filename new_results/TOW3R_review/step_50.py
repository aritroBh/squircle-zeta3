import pickle
from fractions import Fraction as Fr
base = '<project>/research/08_theory/'
H = pickle.load(open(base+'TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def gdiv(x,y):
    den = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den, (x[1]*y[0]-x[0]*y[1])/den)
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def qsi_div(x, y):
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    den = (gmul(C,C)[0]-2*gmul(D,D)[0], gmul(C,C)[1]-2*gmul(D,D)[1])
    numA = (gmul(A,C)[0]-2*gmul(B,D)[0], gmul(A,C)[1]-2*gmul(B,D)[1])
    numB = (gmul(B,C)[0]-gmul(A,D)[0], gmul(B,C)[1]-gmul(A,D)[1])
    return (gdiv(numA, den), gdiv(numB, den))
def qsi_scale(x, c):  # multiply by Gaussian rational c=(re,im)
    A=(x[0],x[2]); B=(x[1],x[3])
    return (gmul(A,c)[0], gmul(B,c)[0], gmul(A,c)[1], gmul(B,c)[1])
for p in [3,4,5]:
    e_pp   = H[(p,p)][0][2*p]
    e_p1p  = H[(p+1,p)][0][2*p+1]
    e_pm1  = H[(p,p-1)][0][2*p-1]
    r1 = qsi_div(e_p1p, e_pp)   # want 2i(4p-1)/(p+1) ?
    r2 = qsi_div(e_p1p, e_pm1)  # want (4p-1)(4p-3)/((p+1)p) ?
    rA = qsi_div(H[(p+1,p-1)][0][2*p], e_pm1)  # want i(4p-3)/(p+1)?
    print(f'p={p}: e00(p+1,p)/e00(p,p) = {r1}   target (0,{2*(4*p-1)}/{p+1})')
    print(f'      e00(p+1,p)/e00(p,p-1) = {r2}  target ({(4*p-1)*(4*p-3)}/{(p+1)*p},0)')
    print(f'      e00(p+1,p-1)/e00(p,p-1) = {rA}  target (0,{4*p-3}/{p+1})')
