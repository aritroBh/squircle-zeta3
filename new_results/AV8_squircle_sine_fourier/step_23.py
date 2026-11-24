src = open('squircle.py').read()
newH = '''
def _gs(y):
    s2 = y*y; ls = log(y); lc2 = log(1-s2); d = s2-1
    g1 = (-2*s2*ls + d*lc2)/(4*d)
    g2 = (2*s2*s2*ls**2 - 2*s2*d*(ls-1)*ls - d*d*lc2)/(4*d*d)
    g3 = (-8*s2**3*ls**3 + 6*s2*s2*d*(2*ls-1)*ls**2 + 6*s2*d*d*(ls-1)*ls
          - 4*s2*d*d*ls**3 + 3*d**3*lc2)/(8*d**3)
    return g1, g2, g3

def H1u(u):
    s = sin(u); c = cos(u)
    g1,g2,g3 = _gs(s)
    H = (1-s*s)**mpf('-0.5')
    return H*g1*c

def H2u(u):
    s = sin(u); c = cos(u)
    g1,g2,g3 = _gs(s)
    H = (1-s*s)**mpf('-0.5')
    return H*(g2+g1*g1)*c

def H3u(u):
    s = sin(u); c = cos(u)
    g1,g2,g3 = _gs(s)
    H = (1-s*s)**mpf('-0.5')
    return H*(g3+3*g1*g2+g1**3)*c
'''
import re
# replace old H defs block (from "def H1u" to "def L1")
start = src.index("def H1u(u):")
end = src.index("def L1(phi):")
src = src[:start] + newH.strip() + "\n\n" + src[end:]
open('squircle.py','w').write(src)
print("patched")
