import pickle, sys
sys.path.insert(0,'.')
import engine14 as E
sigma = pickle.load(open('sigma14.pkl','rb'))
print("sigma(S)_w3 terms (letter idx 0=2,1=alpha,2=u):")
for key in sorted(sigma.keys()):
    print("  ", key, " coef=", E.fev(sigma[key],25), " FT:", sigma[key][:3], "/", sigma[key][3])
