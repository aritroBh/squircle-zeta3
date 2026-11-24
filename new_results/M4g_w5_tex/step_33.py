import json
from fractions import Fraction
import mpmath as mp
mp.dps=100
sq2=mp.sqrt(2)
def li32(x,y,N):
    s,xn,yj,h=mp.mpf(0),x,mp.mpf(1),mp.mpf(0)
    for n in range(2,N+1):
        yj*=y; h+=yj/mp.mpf(n-1)**2; xn*=x
        s+=xn/mp.mpf(n)**3*h
    return s
def li32_slow(x,y,N):
    s=mp.mpf(0)
    for n in range(2,N+1):
        xn=x**n
        t=mp.mpf(0); yj=y
        for m in range(1,n):
            t+=yj/mp.mpf(m)**2; yj*=y
        s+=xn/mp.mpf(n)**3*t
    return s
for x,y in ((sq2-1,1),(1/sq2,1),(1-sq2,1),(sq2-1,-1)):
    a=li32(x,y,3500); b=li32_slow(x,y,500)
    print('li32 x=%s y=%s slow-fast diff %s'%(mp.nstr(x,4),mp.nstr(y,2),mp.nstr(b-a,3)))
LCHI=(mp.zeta(5,mp.mpf(1)/8)+mp.zeta(5,mp.mpf(7)/8)-mp.zeta(5,mp.mpf(3)/8)-mp.zeta(5,mp.mpf(5)/8))/8**5
ref=json.load(open('../M4e_w5_lchi/Lchi8_5.json'))
print(ref.keys() if isinstance(ref,dict) else type(ref))
rs=ref['value'] if isinstance(ref,dict) and 'value' in ref else (ref if isinstance(ref,str) else ref.get('L'))
print('LCHI diff vs certified:', mp.nstr(LCHI-mp.mpf(rs),3))
