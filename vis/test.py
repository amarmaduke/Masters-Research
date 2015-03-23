import numpy as np

def test() :
    test = np.zeros((2, 2, 2))
    test[1, 1, 1] = 1
    return test

t = test()
print(t)