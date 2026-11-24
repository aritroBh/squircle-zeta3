# The 4th identity I guessed may not be a true identity. Check numerically which weight-3 shuffle holds.
import mpmath as mp
mp.mp.dps = 40
def li2v(a,b,x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    tot=mp.mpf(0); inner=mp.mpf(0); xn=x; n=1
    tol=mp.mpf(10)**(-mp.mp.dps+10)
    while True:
        if n>=2:
            m=n-1; inner += y**m/mp.mpf(m)**b
        term = xn*inner/mp.mpf(n)**a
        tot += term
        if n>60 and abs(term) < tol: break
        n+=1; xn*=x
    return tot
x,y = mp.mpf('0.5'), mp.mpf('0.7')
Li1=lambda z: -mp.log(1-z)
Li2=lambda z: mp.polylog(2,z)
Li3=lambda z: mp.polylog(3,z)
cands = {
 'Li12(x,y)+Li21(y,x)+Li3(xy)': li2v(1,2,x,y)+li2v(2,1,y,x)+Li3(x*y),
 'Li21(x,y)+Li12(y,x)+Li3(xy)': li2v(2,1,x,y)+li2v(1,2,y,x)+Li3(x*y),
 'Li12(x,y)+Li21(x,y)': li2v(1,2,x,y)+li2v(2,1,x,y),
}
print("Li2(x)Li1(y) =", mp.nstr(Li2(x)*Li1(y), 30))
for k,v in cands.items():
    print(k, "=", mp.nstr(v,30), " diff:", mp.nstr(v-Li2(x)*Li1(y),6))
# stuffle check for weight 3 convB: Li2(x)Li1(y) =? Li21(x,y)+Li12(y,x)+Li3(xy)? tested above.
