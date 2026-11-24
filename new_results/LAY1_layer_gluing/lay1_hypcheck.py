# LAY1: grid sanity checks of the explicit inequalities used in the hypothesis verifications (not part of the proofs).
import mpmath as mp
mp.mp.dps = 400
g1 = lambda s: mp.log(1+mp.exp(-2*abs(s)))
def L_sq(s, t):   # squircle window exponent, |u|=t|s| <= pi/8
    x = t*abs(s); Lam = mp.atanh(mp.tan(x)); return t*mp.log(1+mp.exp(-2*Lam/t))
def L_sm(s, t):   # soft-max window exponent, |u|=t|s| <= pi/4
    u = t*s
    return t*mp.log(1+mp.exp(-2*mp.sin(abs(u))/t)+mp.exp(-2*mp.cos(u)/t)*(1+mp.exp(-2*abs(mp.sin(u))/t)))
worst_sq = worst_sm = mp.mpf(0)
for t in [mp.mpf(1)/n for n in (2, 5, 10, 50, 200, 1000)]:
    for k in range(1, 400):
        s = (mp.pi/8/t)*k/400
        b = (4*mp.sqrt(2)/3)*t**3*s**3*mp.exp(-2*s)
        worst_sq = max(worst_sq, abs(L_sq(s,t)-t*g1(s))/b)
        s2 = (mp.pi/4/t)*k/400
        b2 = t**3*s2**3*mp.exp(-mp.mpf('1.79')*s2)/3 + 2*t*mp.exp(-mp.sqrt(2)/t)
        worst_sm = max(worst_sm, abs(L_sm(s2,t)-t*g1(s2))/b2)
print("max |L-tau g1| / claimed bound: squircle", mp.nstr(worst_sq,6), " softmax", mp.nstr(worst_sm,6), "(must be <= 1)")
# convexity of the soft-max: min over [0,pi/4] of h+h'' (numerical derivative) for several tau
for t in [mp.mpf(1), mp.mpf('0.5'), mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.05')]:
    h = lambda th: (mp.cos(th)+mp.sin(th))*(1+mp.exp(-2*mp.sin(th)/t)+mp.exp(-2*mp.cos(th)/t)+mp.exp(-2*(mp.cos(th)+mp.sin(th))/t))**t
    mn = min(h(th)+mp.diff(h, th, 2) for th in [mp.pi/4*k/200 for k in range(0, 201)])
    # also check Q >= 0 lower bound formula (sum of nonneg terms) is consistent: report min(h+h'')/(t*h)
    print(f"softmax tau={mp.nstr(t,3)}: min_[0,pi/4] (h+h'') = {mp.nstr(mn,8)}")
# outside-window bound for squircle: |h - m| <= tau * tan(pi/8)^(1/tau)
