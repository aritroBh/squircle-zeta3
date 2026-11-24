import sys, os, pickle
from fractions import Fraction as Fr
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def u(i,j,k,delta):
    return tuple(H[(i,j)][delta][k]) if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta] else (Fr(0),)*4
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gdiv(x,y):
    den = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den, (x[1]*y[0]-x[0]*y[1])/den)
def qsi_div(x, y):
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    C2m2D2 = (gmul(C,C)[0]-2*gmul(D,D)[0], gmul(C,C)[1]-2*gmul(D,D)[1])
    numA = (gmul(A,C)[0]-2*gmul(B,D)[0], gmul(A,C)[1]-2*gmul(B,D)[1])
    numB = (gmul(B,C)[0]-gmul(A,D)[0], gmul(B,C)[1]-gmul(A,D)[1])
    return (gdiv(numA, C2m2D2), gdiv(numB, C2m2D2))
r = qsi_div(u(3,2,5,0), u(2,2,4,0))
print('e00(3,2)/e00(2,2) =', r, '   2i(4p-1)/(p+1) at p=2 =', (Fr(0), Fr(14,3)))
# so rB0d should be?
# r = A + B s2; 2i(4p-1)/(p+1) = (14/3) i -> (A=0,B=(14/3,0))? i = (re0,im1): A=(0, 14/3), B=(0,0).
print('A part:', r[0], 'B part:', r[1])
# check u(2,2,5,1)/e00 = -2-2i as computed
rr = qsi_div(u(2,2,5,1), u(2,2,4,0))
print('data ehat(1,0) ratio:', rr)
