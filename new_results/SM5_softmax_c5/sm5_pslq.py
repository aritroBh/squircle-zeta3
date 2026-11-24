# SM5 item 3: PSLQ of the numerical c5.
import mpmath as mp
mp.mp.dps = 55
c5 = mp.mpf("1.48992283096334707297341790385558242735871286")   # from sm5_richardson_output.txt (about 50 digits trusted)
z=mp.zeta; l2=mp.log(2)
def run(name, basis, labels, dps=48):
    mp.mp.dps = dps
    rel = mp.pslq([c5]+basis, maxcoeff=10**6, maxsteps=10**6)
    print(name, "->", rel, " labels:", ["c5"]+labels)
    return rel
run("narrow", [z(5), z(2)*z(3)], ["z5","z2z3"])
run("narrow+z4 (weight mix check)", [z(5), z(2)*z(3), z(4)], ["z5","z2z3","z4"])
wide = [z(5), z(2)*z(3), z(4)*l2, z(3)*l2**2, z(2)*l2**3, l2**5, mp.polylog(5,mp.mpf(1)/2), mp.polylog(4,mp.mpf(1)/2)*l2]
run("wide weight-5 alternating", wide, ["z5","z2z3","z4 l2","z3 l2^2","z2 l2^3","l2^5","Li5(1/2)","Li4(1/2) l2"], dps=50)
