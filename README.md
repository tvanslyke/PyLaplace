# PyLaplace
Provides a numerical Laplace Transform function to CPython.  Allows for applying the Laplace Transform to arbitrary Python functions.

## Implementation
The implementation uses the (Gauss-Laguerre quadrature)[https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature] to approximate the exact Laplace Transform integral.  Because this method is not always numerically stable, the generated functor is divergent from the true transform on part of the complex plane.

The goal of the implementation was to prioritize the accuracy of the approximation above speed and memory usage.  The roots of the first 100 Laguerre polynomials (needed to compute the numerical approximation) are precomputed in the interest of speed, however.  Higher `order` approximations can be particularly expensive because the implementation of the root-finding algorithm relies on an unproven (at least, I think it's unproven) but useful-in-practice property of Laguerre polynomial roots: that the roots of `n`th Laguerre polynomial interleave the roots of the `n+1`th Laguerre Polynomial.  So in order to compute (without precomputing the roots), say `laplace.LaplaceTransform(math.exp, order = 30)`, the root-finding algorithm would have to recursively bucket-brigade the root computation up from `order=1`.  This is why theres a `precomputed_laguerre_roots.cpp`  with a LOT of numbers in it. :)


## Example
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
```
