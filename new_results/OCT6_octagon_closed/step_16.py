import pickle, sympy as sp
simple = pickle.load(open('mons_simple.pkl','rb'))
r2 = sp.sqrt(2)
for key, terms in simple:
    if key in [(0,0,1,2,0,0),(0,0,1,0,0,2)]:
        for n,q in terms:
            if n==2:
                print(key, 'raw q:', q)
                # force into p + q*r2 form
                qq = sp.expand(q)
                print('collect:', sp.collect(qq, r2))
                print('numeric:', sp.N(qq, 30))
