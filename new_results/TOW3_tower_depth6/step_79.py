import sys, os, pickle
from fractions import Fraction as Fr
from math import comb
base = '<project>/research/08_theory'
H = pickle.load(open(f'{base}/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def u(i,j,k,delta):
    return tuple(H[(i,j)][delta][k]) if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta] else (Fr(0),)*4
def qsi_add(x,y): return tuple(a+b for a,b in zip(x,y))
def qsi_scale(x,f): return tuple(a*f for a in x)
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
def e_ratio(key, delta, s):
    i,j = key
    k0 = (2*i+1) if i==j else i+j
    out = (Fr(0),)*4
    for t in range(0, s+1):
        out = qsi_add(out, qsi_scale(u(i,j,k0+t,delta), Fr((-1)**(s-t)*comb(s,t))))
    return qsi_div(out, u(i,j,k0,0))
parent=(3,0); d=2; n=3
print('e_ratio(2,1):', e_ratio(parent,2,1))
print('e_ratio(2,0):', e_ratio(parent,2,0))
print('e_ratio(1,0):', e_ratio(parent,1,0))
print('e00:', u(3,0,3,0))
# expected beta_{2,2} = -3 * 2Re[E * (-1)*e_ratio(2,1) * w_1], w_1 = s2(-1/2-i/2)
E = ((u(3,0,3,0)[0], u(3,0,3,0)[2]), (u(3,0,3,0)[1], u(3,0,3,0)[3]))
print('E =', E)
erg = e_ratio(parent,2,1)[0]
print('erg =', erg)
# Y = -erg*w1: w1 rat part (0,0), s2 part (-1/2,-1/2)
# E = (0 + (1/4) s2) i.e. Er=(0,0), Es=(1/4,0)
# E*Y = (s2/4) * [-erg * s2 * (-1/2-i/2)] = (2/4)(-erg)(-1/2-i/2) = erg(1/4+... wait
# -erg*w1 = -(-1/3)*s2*(-1/2-i/2) = (1/3)(-1/2-i/2) s2 = (-1/6-i/6) s2
# E*Y = (s2/4)*(-1/6-i/6)*s2 = (1/2)(-1/6-i/6) = -1/12 - i/12
# 2Re = -1/6; beta = -3 * (-1/6)?? = 1/2... 
# Hmm that gives 1/2! So beta = -(i+1)*2Re[E Y] with 2Re[E Y] = -1/6 gives 1/2.
# But note: 2Re[E Y] = -1/6 means beta = -3 * (-1/6) = 1/2 ✓
print('expected beta = 1/2')
