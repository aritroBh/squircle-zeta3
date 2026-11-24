import os
print("CWD:", os.getcwd())
print("CWD exists:", os.path.exists(os.getcwd()))
# test write
try:
    open('_test_write.tmp','w').write('x')
    print("write OK"); os.remove('_test_write.tmp')
except Exception as e:
    print("write FAIL:", e)
# check the a1=1.00441 value in extra
import pickle
extra=pickle.load(open('gamma_extra.pkl','rb'))
for r in extra:
    print("extra a1:", repr(float(r[0])))
