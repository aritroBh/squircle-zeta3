import mpmath as mp
mp.mp.dps=80
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a)
def g1(x): return mp.log((x-b)*(x+a)/(2*x))
def g2(x): return mp.log((x*x+a*a)*(x*x+b*b)/((x-b)**2*(x+a)**2))
def tg3(x): return mp.log((x*x+b*b)/(x*x+a*a))
def dg1(x): return 1/(x-b)+1/(x+a)-1/x
def dg2(x): return 2*x/(x*x+a*a)+2*x/(x*x+b*b)-2/(x-b)-2/(x+a)
def dtg3(x): return 2*x/(x*x+b*b)-2*x/(x*x+a*a)
def OmA(x): return 1/x - r2/2/(x-b) + r2/2/(x+a)
def OmB(x): return -1/x + 1/(x-b) + 1/(x+a)
def Q(f): return mp.quad(f, [1, 8, 64, mp.inf])
G1X=Q(lambda x: OmA(x)*g1(x)*g2(x))
G2X=Q(lambda x: (dg1(x)*g2(x)+g1(x)*dg2(x))/(x-b))
G3X=Q(lambda x: (dg1(x)*g2(x)+g1(x)*dg2(x))/(x+a))
H1X=Q(lambda x: OmB(x)*g1(x)*tg3(x))
H2X=Q(lambda x: (dg1(x)*tg3(x)+g1(x)*dtg3(x))/(x-b))
H3X=Q(lambda x: (dg1(x)*tg3(x)+g1(x)*dtg3(x))/(x+a))
H4X=Q(lambda x: g1(x)/(x-b)**2)
H5X=Q(lambda x: g1(x)/(x+a)**2)
V1X = -G1X - b*G2X + a*G3X
V3X = -H1X + (b*r2/2)*H2X + (a*r2/2)*H3X + r2*b*L*H4X + r2*a*L*H5X
print("V1X:", mp.nstr(V1X,25), " V3'X:", mp.nstr(V3X,25))

# y-form pieces
def hg1(y): return mp.log((1-b*y)*(1+a*y)/(2*y))
def hg2(y): return mp.log((1+a*a*y*y)*(1+b*b*y*y)/((1-b*y)**2*(1+a*y)**2))
def hg3(y): return mp.log((1+b*b*y*y)/(1+a*a*y*y))
def Qy(f): return mp.quad(f, [0, mp.mpf('0.3'), mp.mpf('0.7'), 1])
G1y=Qy(lambda y: (1/y - r2*b/2/(1-b*y) - r2*a/2/(1+a*y))*hg1(y)*hg2(y))
G2y=Qy(lambda y: hg1(y)*hg2(y)/(1-b*y)**2)
G3y=Qy(lambda y: hg1(y)*hg2(y)/(1+a*y)**2)
H1y=Qy(lambda y: (1/y + b/(1-b*y) - a/(1+a*y))*hg1(y)*hg3(y))
H2y=Qy(lambda y: hg1(y)*hg3(y)/(1-b*y)**2)
H3y=Qy(lambda y: hg1(y)*hg3(y)/(1+a*y)**2)
H4y=Qy(lambda y: hg1(y)/(1-b*y)**2)
H5y=Qy(lambda y: hg1(y)/(1+a*y)**2)
for nm,X,Y in [("G1",G1X,G1y),("G2",G2X,G2y),("G3",G3X,G3y),("H1",H1X,H1y),("H2",H2X,H2y),("H3",H3X,H3y),("H4",H4X,H4y),("H5",H4X,H4y)]:
    print(nm, "X:", mp.nstr(X,15), " y:", mp.nstr(Y,15), " diff:", mp.nstr(abs(X-Y),4))
