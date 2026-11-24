exec(open('common.py').read())
C4, Q7 = get_C4_Q7()
rev = [(a1,a5),(a2,a4),(a4,a2),(a5,a1)]  # a3 fixed
print("C4 rev-invariant:", sp.expand(C4 - C4.subs(rev, simultaneous=True))==0)
print("Q7' rev-invariant:", sp.expand(Q7 - Q7.subs(rev, simultaneous=True))==0)
# Also as dicts
def symrev(p):
    return sp.expand(p.subs(rev, simultaneous=True))
# earlier failure: subs(dict) maybe did sequential? test
d1 = C4.subs({a1:a5,a2:a4,a4:a2,a5:a1})
print("dict subs diff:", sp.expand(C4-d1)==0)
