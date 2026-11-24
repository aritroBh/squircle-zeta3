from cor8r_num import *
setdps(50)
x, y = mp.mpf('0.3'), mp.mpf('0.4')

def G_0a0b_series(a, b, K=200):
    """G(0,a,0,b;1) = sum_{k2>k1>=1} a^{-(k2-k1)} b^{-k1}/(k1^2 k2^2)   [derived by expanding the integral]"""
    a = mp.mpf(a); b = mp.mpf(b)
    s = mp.mpf(0)
    for k2 in range(2,K):
        for k1 in range(1,k2):
            s += a**(-(k2-k1)) * b**(-k1)/(k1**2 * k2**2)
    return s

def G_a00b_series(a, b, K=200):
    """G(a,0,0,b;1) = sum_{k2>k1>=1} a^{-(k2-k1)} b^{-k1}/(k1^3 k2)   [derived by expanding the integral]"""
    a = mp.mpf(a); b = mp.mpf(b)
    s = mp.mpf(0)
    for k2 in range(2,K):
        for k1 in range(1,k2):
            s += a**(-(k2-k1)) * b**(-k1)/(k1**3 * k2)
    return s

L22 = Li_ab(2,2,x,y)
L31 = Li_ab(3,1,x,y)
G22 = G_0a0b_series(1/y, 1/(x*y))
G31 = G_a00b_series(1/y, 1/(x*y))
print("Li_{2,2}(x,y) series      :", mp.nstr(L22,40))
print("G(0,1/y,0,1/(xy);1) series:", mp.nstr(G22,40))
print("residual:", mp.nstr(abs(L22-G22),5))
print()
print("Li_{3,1}(x,y) series      :", mp.nstr(L31,40))
print("G(1/y,0,0,1/(xy);1) series:", mp.nstr(G31,40))
print("residual:", mp.nstr(abs(L31-G31),5))
# convB check: convB-Li_{a,b}(x,y) = G(0^{a-1},1/x,0^{b-1},1/(xy);1) should equal convA-Li_{b,a}(y,x)
G22b = G_0a0b_series(1/x, 1/(x*y))
print()
print("G(0,1/x,0,1/(xy);1) series:", mp.nstr(G22b,40))
print("convA Li_{2,2}(y,x)        :", mp.nstr(Li_ab(2,2,y,x),40))
print("residual:", mp.nstr(abs(G22b-Li_ab(2,2,y,x)),5))
