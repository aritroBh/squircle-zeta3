import re, sys, pickle
sys.path.insert(0,'lib')
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application
import tow2lib as T2, tow4lib as T4
d,m=T2.d,T2.m
loc={'d':d,'m':m,'I':sp.I}
cache={k:sp.sympify(v,locals=loc) for k,v in pickle.load(open('lib/forced_F78.pkl','rb')).items()}
md=open('../../proof_depth8.md').read()
tr=standard_transformations+(implicit_multiplication_application,)
def P(txt, pil):
    t=txt.replace('^','**')
    t=re.sub(r'\bPi_(\d)', lambda mo: f'({pil[int(mo.group(1))]})', t)
    t=re.sub(r'(?<![A-Za-z])i(?![A-Za-z])','I',t)
    return parse_expr(t, local_dict={'d':d,'m':m,'I':sp.I}, transformations=tr)
pil={3:'(2*d-7)*(2*d-5)*(2*d-3)',4:'(2*d-9)*(2*d-7)*(2*d-5)*(2*d-3)'}
# grab the block between "The new closed forms" and "[The tops"
blk=md[md.index('The new closed forms'):md.index('[The tops')]
blk=blk.replace('\n',' ')
items=re.findall(r'F\((\d),(\d)\)\s*=\s*(.+?)(?=\s{2}F\(\d,\d\)\s*=|\s*$)', blk)
print('parsed displayed forms:', [(a,b) for a,b,_ in items])
for a,b,txt in items:
    key=(int(a),int(b))
    txt=txt.strip().rstrip('.').strip()
    # drop trailing commentary in parentheses after the formula, e.g. "(zigzag top, A_4 = 5)"
    txt=re.sub(r'\((?:zigzag|[^()]*A_\d)[^()]*\)\s*$','',txt).strip()
    txt=re.sub(r'\s*\(?-?1\)?\^4 A_5/Pi_4 with A_5 = 16\)?','',txt)
    try:
        e=P(txt,pil)
    except Exception as ex:
        print(key,'PARSE FAIL',ex, txt[:80]); continue
    diff=sp.cancel(sp.together(e-cache[key]))
    print(key,'displayed == cache:', diff==0, ('' if diff==0 else '   diff = '+str(sp.factor(diff))[:100]))
# a_7 displayed numerator
m7=re.search(r'N_7 =\s*(.+?)\n\s*a_8',md,re.S)
n7=m7.group(1).replace('\n',' ')
e7=P(n7,pil)/(40320*(2*d-7)*(2*d-5)*(2*d-3))
print('a_7 displayed == cache F(7,0):', sp.cancel(sp.together(e7-cache[(7,0)]))==0)
