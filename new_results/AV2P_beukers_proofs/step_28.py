from mpmath import mp, mpf, log, pi, nprint, quad
mp.dps = 50
x = mpf(1)/16
def A(x):
    return mp.asin(2*mp.sqrt(x))/(2*mp.sqrt(x))
# B-C = int [A(x t^2) - A(x t)]/(1-t) dt
BC = quad(lambda t: (A(x*t*t)-A(x*t))/(1-t), [0, 1])
print("int [A(xt^2)-A(xt)]/(1-t) =", nprint(BC, 40))
print("B-C (series) =", nprint(mpf('0.0502870873772545765534176332325716727641098374883945464791628')-mpf('0.0743778361833922484294700798458102948671232805137595100031871'), 40))
W = 2*quad(lambda t: mp.asin(t/2)/(t*(1+t)), [0,1])
print("W =", nprint(W, 40))
# pieces
P1 = quad(lambda t: 2*mp.asin(t/2)/t - 2*mp.asin(mp.sqrt(t)/2)/mp.sqrt(t), [0,1])
print("int [2asin(t/2)/t - 2asin(sqrt t/2)/sqrt t]/(1-t)... full:", nprint(quad(lambda t: (2*mp.asin(t/2)/t - 2*mp.asin(mp.sqrt(t)/2)/mp.sqrt(t))/(1-t), [0,1]), 40))
# where is the discrepancy? compute int A(xt^2)/(1-t) - A(xt)/(1-t) split:
J1 = quad(lambda t: (2*mp.asin(t/2)/t - pi/3)/(1-t), [0,1])  # subtract limit for convergence
J2 = quad(lambda t: (2*mp.asin(mp.sqrt(t)/2)/mp.sqrt(t) - pi/3)/(1-t), [0,1])
print("J1 - J2 =", nprint(J1-J2, 40))
