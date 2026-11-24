import sys, pickle
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import tow4lib as T4
work = '<project>/research/08_theory/TOW4_tower_depth8'
a7 = sp.sympify(pickle.load(open(work+'/a7_closed.pkl','rb')), locals={'d':T4.d,'m':T4.m,'I':sp.I})
a8 = sp.sympify(pickle.load(open(work+'/a8_closed.pkl','rb')), locals={'d':T4.d,'m':T4.m,'I':sp.I})
print("arho(7)==a7_closed:", sp.cancel(T4.arho(7)-a7)==0)
print("arho(8)==a8_closed:", sp.cancel(T4.arho(8)-a8)==0)
# update forced78.py to use tow4lib tau
src = open(work + '/forced78.py').read()
src = src.replace("import tow2lib as T2\nd, m = T2.d, T2.m\nI = sp.I\ntau1 = T2.tau1f; tau2 = T2.tau2f",
                  "import tow4lib as T4\nimport tow2lib as T2\nd, m = T2.d, T2.m\nI = sp.I\ntau1 = T4.tau1f; tau2 = T4.tau2f")
open(work + '/forced78.py','w').write(src)
print("forced78.py updated to tow4lib tau")
