
# step_01.py -- fast convB MPL series + G-word series formulas; verify word identification
import mpmath as mp
mp.mp.dps = 60

def li_convB(a, b, x, y):
    # Sum_{n>m>=1} x^n y^m / (n^a m^b); converges for |x|<1, |x*y|<=1 (x*y<1 safe)
    x = mp.mpf(x); y = mp.mpf(y)
    total = mp.mpf(0); m = 0
    partial = mp.mpf(0); xp = mp.mpf(1)
    Li_a_x = mp.polylog(a, x)
    tol = mp.mpf(10)**(-mp.mp.dps+8)
    while True:
        m += 1
        xp *= x
        partial += xp / mp.power(m, a)
        tail = Li_a_x - partial
        term = mp.power(y, m) / mp.power(m, b) * tail
        total += term
        if m > 40 and abs(term) < tol:
            break
    return total

def G_word_series(m1, c1, m2, c2, m3=None, c3=None):
    # G_phys(0^{m1},c1,0^{m2},c2;1) = sum_{k>j>=1} c1^{-(k-j)} c2^{-j} / (j^{m2+1} k^{m1+1})
    # depth-3: G_phys(0^{m1},c1,0^{m2},c2,0^{m3},c3;1)
    #        = sum_{k>j>i>=1} c1^{-(k-j)} c2^{-(j-i)} c3^{-i} / (i^{m3+1} j^{m2+1} k^{m1+1})
    tol = mp.mpf(10)**(-mp.mp.dps+8)
    if m3 is None:
        total = mp.mpf(0); j = 0
        c2p = mp.mpf(1)
        while True:
            j += 1; c2p /= c2
            # inner: sum_{k>j} c1^{-(k-j)}/k^{m1+1} = c1^j * (Li_{m1+1}(1/c1) - sum_{k<=j} c1^{-k}/k^{m1+1})
            partial = mp.mpf(0); c1p = mp.mpf(1)
            for k in range(1, j+1):
                c1p /= c1
                partial += c1p / mp.power(k, m1+1)
            tailK = mp.polylog(m1+1, 1/c1) - partial
            term = c2p/mp.power(j, m2+1) * (c1**j) * tailK
            total += term
            if j > 40 and abs(term) < tol:
                break
        return total
    else:
        total = mp.mpf(0); i = 0; c3p = mp.mpf(1)
        while True:
            i += 1; c3p /= c3
            # sum over j>i of c2^{-(j-i)}/j^{m2+1} * [tail over k]
            inner = mp.mpf(0); j = i; c2p = c2**i
            while True:
                j += 1; c2p /= c2
                partial = mp.mpf(0); c1p = mp.mpf(1)
                for k in range(1, j+1):
                    c1p /= c1
                    partial += c1p / mp.power(k, m1+1)
                tailK = mp.polylog(m1+1, 1/c1) - partial
                termj = c2p/mp.power(j, m2+1) * (c1**j) * tailK
                inner += termj
                if j > i+40 and abs(termj) < tol:
                    break
            term = c3p/mp.power(i, m3+1) * inner
            total += term
            if i > 40 and abs(term) < tol:
                break
        return total

x, y = mp.mpf('0.5'), mp.mpf('0.7')
print("point x=0.5 y=0.7")
print("Li22 convB :", mp.nstr(li_convB(2,2,x,y), 40))
print("G(0,1/x,0,1/xy):", mp.nstr(G_word_series(1,1/x,1,1/(x*y)), 40))
print("Li31 convB :", mp.nstr(li_convB(3,1,x,y), 40))
print("G(0,0,1/x,1/xy):", mp.nstr(G_word_series(2,1/x,0,1/(x*y)), 40))
print("Li13 convB :", mp.nstr(li_convB(1,3,x,y), 40))
print("G(1/x,0,0,1/xy):", mp.nstr(G_word_series(0,1/x,2,1/(x*y)), 40))
print("Li21 convB :", mp.nstr(li_convB(2,1,x,y), 40))
print("G(0,1/x,1/xy)  :", mp.nstr(G_word_series(1,1/x,0,1/(x*y)), 40))
print("Li12 convB :", mp.nstr(li_convB(1,2,x,y), 40))
print("G(1/x,0,1/xy)  :", mp.nstr(G_word_series(0,1/x,1,1/(x*y)), 40))
print("Li11 convB :", mp.nstr(li_convB(1,1,x,y), 40))
print("G(1/x,1/xy)    :", mp.nstr(G_word_series(0,1/x,0,1/(x*y)), 40))

# shuffle check: Li22(x,y)+Li22(y,x)+Li4(xy) = Li2(x)Li2(y)
lhs = li_convB(2,2,x,y) + li_convB(2,2,y,x) + mp.polylog(4, x*y)
rhs = mp.polylog(2,x)*mp.polylog(2,y)
print("shuffle22 residual:", mp.nstr(lhs-rhs, 8))
# stuffle: Li3(x)Li1(y) = Li31(x,y) + Li13(y,x) + Li4(xy)
lhs2 = li_convB(3,1,x,y) + li_convB(1,3,y,x) + mp.polylog(4,x*y)
rhs2 = mp.polylog(3,x)*(-mp.log(1-y))
print("stuffle31 residual:", mp.nstr(lhs2-rhs2, 8))
