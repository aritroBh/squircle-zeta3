import os
src = open('dilog_exact.py').read()
i = src.find('def float_exact')
print(src[i:i+400])
import sys
sys.dont_write_bytecode = True
sys.path.insert(0,'.')
import mpmath as mp, pickle
mp.mp.dps = 120
tails = pickle.load(open('m4_tails_exact.pkl','rb'))
W3exact = tails['W3p'].num()
print(mp.nstr(W3exact, 110))
