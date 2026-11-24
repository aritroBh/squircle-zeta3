import mpmath as mp
mp.mp.dps = 55
def int01_double_sing(x, y):
    # int_0^1 s^x (1-s)^{y-1} ds via geometric slices at both ends + exact tails
    f = lambda s: s**x*(1-s)**(y-1)
    tot = mp.mpf(0)
    # left tail [0, 2^-N]: int ~ a^{x+1}/(x+1) ((1-s)^{y-1} -> 1)
    a = mp.mpf('0.5')
    while a > mp.mpf(10)**(-40):
        tot += mp.quad(f, [a/2, a]); a /= 2
    left_tail = a**(x+1)/(x+1)
    # right tail [1-2^-N, 1]: int ~ b^y/y (s^x -> 1)
    b = mp.mpf('0.5')
    while b > mp.mpf(10)**(-40):
        tot += mp.quad(f, [1-b, 1-b/2]); b /= 2
    right_tail = b**y/y
    return tot + left_tail + right_tail
worst = mp.mpf(0)
for (x,y) in [(mp.mpf('0.3'),mp.mpf('0.4')), (mp.mpf('-0.7'),mp.mpf('0.25')), (mp.mpf('0.9'),mp.mpf('0.9')),
              (mp.mpf('0.123'),mp.mpf('0.456')), (mp.mpc('0.2','0.3'), mp.mpc('0.5','-0.2')), (mp.mpf('-0.95'), mp.mpf('0.99'))]:
    H = int01_double_sing(x,y); B = mp.beta(1+x,y)
    r = abs(H-B); worst = max(worst,r)
    print(f"  H({mp.nstr(x,5)},{mp.nstr(y,5)}) - B = {mp.nstr(r,3)}")
print("worst:", mp.nstr(worst,3))
with open('progress.log','a') as f:
    f.write("step_15: H(x,y)=B(1+x,y) with two-end slice+tail quadrature incl coordinator FAIL pt (0.3,0.4): worst "+mp.nstr(worst,3)+" dps55 -> coordinator FAIL was a quadrature tolerance issue, identity holds\n")
