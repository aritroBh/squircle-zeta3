import mpmath as mp
mp.dps=450
sq2=mp.sqrt(2)
EPS=mp.mpf(10)**(-440)
def li32P(x,y):
    tot,H,xn,ym,n=mp.mpf(0),mp.mpf(0),mp.mpf(1),mp.mpf(1),1
    while True:
        xn*=x; tot+=xn*H/mp.mpf(n)**3
        if n>50 and abs(xn)<EPS: return tot
        ym*=y; H+=ym/mp.mpf(n)**2; n+=1
def li32M(x,y):
    N=int(440/float(-mp.log10(abs(x))))+80
    s,xn,yj,h=mp.mpf(0),x,mp.mpf(1),mp.mpf(0)
    for n in range(2,N+1):
        yj*=y; h+=yj/mp.mpf(n-1)**2; xn*=x
        s+=xn/mp.mpf(n)**3*h
    return s
b=sq2-1
for x,y in ((b,1),(1/sq2,1),(-b,1),(b,-1)):
    p=li32P(x,y); m=li32M(x,y)
    print('x=%s y=%2s P=%s'%(mp.nstr(x,4),mp.nstr(y,2),mp.nstr(p,20)))
    print('             M=%s  P-M=%s'%(mp.nstr(m,20),mp.nstr(p-m,3)))
