
import pickle, os
import sympy as sp
O5 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'oct5')
def load(name):
    with open(os.path.join(O5, name), 'rb') as f:
        return pickle.load(f)
