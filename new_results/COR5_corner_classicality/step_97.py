import mpmath as mp
mp.mp.dps = 200
def Li(s, x): return mp.polylog(s, x)
v = 6*Li(3, mp.mpf(-2)) + 2*Li(3, mp.mpf(1)/2) - Li(3, mp.mpf(-1)/2) - mp.mpf(3)/4*mp.zeta(3)
print("6Li3(-2)+2Li3(1/2)-Li3(-1/2)-(3/4)zeta3 =", mp.nstr(v, 30))
print("|v| =", mp.nstr(abs(v), 3))
# Also the (2,2) cobracket: identity's delta_22 must vanish (classical Li4 has delta_22=0)
# delta_22(H): symbol (e1,e2,e1,e1)-(e1^4): wedge(e1,e2)^wedge(e1,e1)=0; so delta_22=0. RHS: 0. OK.
