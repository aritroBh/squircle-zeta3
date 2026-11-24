import sys, os, pickle
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
H = pickle.load(open(f'{base}/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def u(i,j,k,delta):
    return tuple(H[(i,j)][delta][k]) if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta] else (Fr(0),)*4
print('e00_{3,0} = u(3,0,3,0) =', u(3,0,3,0))
print('e00_{4,0} = u(4,0,4,0) =', u(4,0,4,0))
print('e00_{2,0} = u(2,0,2,0) =', u(2,0,2,0))
# u(3,0,3,2) (depth 2, k=3) and u(3,0,4,2):
print('h(2) k=3:', u(3,0,3,2), ' k=4:', u(3,0,4,2))
# e(2,1) ratio = u(4,2)-u(3,2) over e00
def qsi_sub(x,y): return tuple(a-b for a,b in zip(x,y))
print('diff h(2):', qsi_sub(u(3,0,4,2), u(3,0,3,2)))
