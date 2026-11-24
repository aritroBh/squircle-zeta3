# step_75: test Kummer-Spence 8-term relation for the trilogarithm on L3
import mpmath as mp
mp.mp.dps = 80
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(abs(1-x))
def KS(x, y):
    x, y = mp.mpf(x), mp.mpf(y)
    return (L3(x) - L3(y) + L3(y/x) + L3((1-x)/(1-y)) - L3((1-x)/(1-x*y)) 
            + L3((1-y)/(1-x*y)) - L3(x*(1-y)/(1-x*y)) + L3(y*(1-x)/(1-x*y)))
for x,y in [(mp.mpf('0.3'),mp.mpf('0.6')), (mp.mpf('0.2'),mp.mpf('0.5')), (mp.mpf('0.4'),mp.mpf('0.7')), ('2','3'), ('-2','0.3')]:
    try:
        print(f"KS({x},{y}) =", mp.nstr(KS(x,y), 10))
    except Exception as e:
        print(x, y, "ERR", e)
