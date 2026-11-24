from towlib import *
x = (IQ-BQ)**(-2) + (IQ+AQ)**(-2) - IQ**(-2)
print('raw sum:', x)
print('expected DL1\'(i) = (-1-i)/2 ->', QSI(Fr(-1,2),0,Fr(-1,2),0))
# check each term
print((IQ-BQ)**(-2), (IQ+AQ)**(-2), IQ**(-2))
