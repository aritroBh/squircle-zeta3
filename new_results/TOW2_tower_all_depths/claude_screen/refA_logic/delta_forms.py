import sys; sys.sys_argv=sys.argv; sys.argv=['x','0']
sys.path.insert(0,'.')
import sympy as sp
src = open('prove_E3b34.py').read().split('import sys\ndelta = int')[0]
ns = {}
exec(src, ns)
p = ns['p']; Delta = ns['Delta']
D3 = sp.factor(sp.cancel(Delta(3))); D4 = sp.factor(sp.cancel(Delta(4)))
print('Delta(3) =', D3)
print('Delta(4) =', D4)
doc3 = -(2*p-3)*(32*p**3-96*p**2+130*p-69)/(12*(4*p-3))
doc4 = sp.I*(512*p**6-3840*p**5+12992*p**4-25600*p**3+30722*p**2-20789*p+6054)/(48*(4*p-5)*(4*p-3))
print('doc Delta(3) matches:', sp.cancel(doc3 - D3) == 0)
print('doc Delta(4) matches:', sp.cancel(doc4 - D4) == 0)
print('Delta(0..2):', [sp.factor(sp.cancel(Delta(k))) for k in range(3)])
