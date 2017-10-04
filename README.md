# PyLaplace
Provides a numerical Laplace Transform function to CPython.  Allows for applying the Laplace Transform to arbitrary Python functions.

```
>>> from laplace import LaplaceTransform as LT
>>> from math import sin
>>> exact = lambda s: 1.0 / (1.0 + s ** 2)
>>> F = LT(sin, order = 30)
>>> for s in range(2, 10):
...     print("exact = %f, approx = %f" % (exact(s).real, F(s).real)) 
... 
exact = 0.200000, approx = 0.250000
exact = 0.100000, approx = 0.111111
exact = 0.058824, approx = 0.062500
exact = 0.038462, approx = 0.040000
exact = 0.027027, approx = 0.027778
exact = 0.020000, approx = 0.020408
exact = 0.015385, approx = 0.015626
exact = 0.012195, approx = 0.012348

