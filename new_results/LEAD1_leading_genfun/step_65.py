import mpmath, sympy as sp
mpmath.mp.dps = 80
from mpmath import mpf
keysum_total = mpf('3969.902649561077265952619857564563777846')
oddint = mpf('181527.0001855750653060593112289602953257')
X = sp.symbols('X'); ser = sp.series((1+sp.sinh(X))/sp.cosh(X),X,0,10).removeO()
delta9 = mpf(str(sp.Rational(ser.coeff(X,9))))
print('keysum/delta_9 =', mpmath.nstr(keysum_total/delta9, 40))
print('oddint         =', mpmath.nstr(oddint, 40))
print('diff           =', mpmath.nstr(keysum_total/delta9 - oddint, 6))
with open('progress.log','a') as f:
    f.write('step_61: corrected: odd-d collapse at d=9 with real kernel data confirmed: keysum/delta_9 == oddint to printed diff\n')
