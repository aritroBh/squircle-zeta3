import sys, math, time, os
sys.dont_write_bytecode=True
import flint
import mpmath as mp
mp.mp.dps = 800  # correct way
# planted relation control: 1, sqrt2, 1+sqrt2 at 500 digits in a small lattice
t0=time.monotonic()
n=3
C=mp.mpf(10)**500
plant=[mp.mpf(1), mp.sqrt(2), 1+mp.sqrt(2)]
rows=[[int(i==j) for j in range(n)]+[int(mp.nint(C*v))] for i,v in enumerate(plant)]
red=flint.fmpz_mat(rows).lll(gram='exact')
for i in range(n):
    print([int(red[i,j]) for j in range(n+1)])
print('planted control OK, %.2fs'%(time.monotonic()-t0))
# check w4_symbol Q and M5a n6_values
sys.path.insert(0,'.')
from w4_symbol import Q,R,A,B
print('R=',R,'A=',A,'B=',B)
print('Q(1/2, -3/8)=',Q(1,2) if False else Q(__import__('fractions').Fraction(1,2), __import__('fractions').Fraction(-3,8)))
p='<project>/research/08_theory/M5a_next_coefficient/n6_values.json'
print('M5a n6_values exists:', os.path.exists(p))
import json
print(json.load(open(p)))
