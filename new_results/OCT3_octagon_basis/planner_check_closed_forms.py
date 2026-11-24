import json, mpmath as mp
mp.mp.dps=200
v=[mp.mpf(x) for x in json.load(open('D15b_polygon_law/vals_m2_8_200.json'))['4']]
s2=mp.sqrt(2); u=mp.sqrt(2+s2); L=mp.log(2); M=mp.log(1+s2); pi=mp.pi
P0=8*(s2-1)*mp.sqrt(2-s2)
P1=u*((40-28*s2)*L+(24-16*s2)*M)
P2=u*((-mp.mpf(7)/30+13*s2/60)*pi**2+(123-87*s2)*L**2+(88-62*s2)*L*M+(mp.mpf(104)/5-68*s2/5)*M**2
      +(-mp.mpf(4)/5+3*s2/5)*mp.polylog(2,12*s2-16)+(-80+56*s2)*L+(-48+32*s2)*M)
for k,p in enumerate([P0,P1,P2]): print(f'P^({k})(1): file - formula = {mp.nstr(v[k]-p,5)}')
C4=192*mp.sin(pi/8)**4*mp.cos(pi/8)*(3*mp.sin(pi/8)**2-1); print('C4 - (54-39√2)u =',mp.nstr(C4-(54-39*s2)*u,5))
