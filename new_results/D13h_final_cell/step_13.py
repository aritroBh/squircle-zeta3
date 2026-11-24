src = open('step_10.py').read()
old = """lamRreg = (pi*lo(R)**2/2)*g1R/(2*c2reg_up).sqrt() - a*up(Y**arb("1.5"))*Q1abs/w_J
lamLreg = (pi*lo(R)**2/2)*g1L/(2*c2reg_up).sqrt() - a*up(Y**arb("1.5"))*Q1abs/w_J"""
new = """yJ32 = y_J**arb("1.5")
lamRreg = (pi*lo(R)**2/2)*g1R/(2*c2reg_up).sqrt() - a*yJ32*Q1abs/w_J
lamLreg = (pi*lo(R)**2/2)*g1L/(2*c2reg_up).sqrt() - a*yJ32*Q1abs/w_J"""
assert old in src
src = src.replace(old, new)
# also remove the junk double-assignment of g1R
old2 = """g1R = lo(rU/rt2)*lo(mGhp_num(ERa, T))/up(Gh_g(ERa, T))**2
g1R = lo(iv(1,1+ (delta_max/2)))*0 + (1/rt2)*lo(mGhp_num(ERa, T))/up(Gh_g(ERa, T))**2 * lo(R)
g1L = (1/rt2)*lo(mGhp_num(ELa, T))/up(Gh_g(ELa, T))**2 * lo(R)"""
new2 = """g1R = (1/rt2)*lo(mGhp_num(ERa, T))/up(Gh_g(ERa, T))**2 * lo(R)
g1L = (1/rt2)*lo(mGhp_num(ELa, T))/up(Gh_g(ELa, T))**2 * lo(R)"""
assert old2 in src
src = src.replace(old2, new2)
open('step_10.py','w').write(src)
import subprocess, sys
r = subprocess.run([sys.executable, 'step_10.py'], capture_output=True, text=True, timeout=300)
print(r.stdout)
print("STDERR:", r.stderr[-3000:] if r.stderr else "(none)")
