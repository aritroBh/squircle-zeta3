import re, sys
sys.path.insert(0, '.')
import sympy as sp
import check_data as CD
import tow3lib as T3
txt = open(CD.MD).read()
d, m = T3.d, T3.m
def Pi_(n): return T3.Pi(n)
def grab(pat):
    mt = re.search(pat, txt, re.S)
    assert mt, pat
    return CD.parse_poly(mt.group(1))
forms = {
 (2,1): -1/Pi_(1),
 (3,1): grab(r'e\(3,1\) = (i \(d-3\))/Pi_1')/Pi_(1),
 (4,1): grab(r'e\(4,1\) = \((.*?)\)/Pi_2')/Pi_(2),
 (4,2): sp.Integer(2)/Pi_(2),
 (5,1): grab(r'e\(5,1\) = \((.*?)\)/\(12 Pi_2\)')/(12*Pi_(2)),
 (5,2): grab(r'e\(5,2\) = (-2 i \(d-4\))/Pi_2')/Pi_(2),
 (6,1): grab(r'e\(6,1\) = \((.*?)\)/\(96 Pi_3\)')/(96*Pi_(3)),
 (6,2): grab(r'e\(6,2\) = \((.*?)\)/\(2 Pi_3\)')/(2*Pi_(3)),
 (6,3): sp.Integer(-5)/Pi_(3),
}
for (dl, s), f in forms.items():
    diff = sp.cancel(sp.together(f - T3.e(dl, s)))
    print(f'md e({dl},{s}) - tow3lib.e({dl},{s}) == 0 : {diff == 0}')
# also zigzag A_n and tau table sanity
import tow2lib as T
print('zigzag A_0..A_5:', [T.A[n] for n in range(6)])
print('check_tau (n<=7):', T.check_tau())
