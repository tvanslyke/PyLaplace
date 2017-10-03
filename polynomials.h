#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include <boost/math/tools/polynomial.hpp>
#include "big_numbers.h"
template <class T>
using poly_t = boost::math::tools::polynomial<T>;


template <class T>
auto evaluate_polynomial(const poly_t<T> & p, const T & arg)
{
	T argpow = 1;
	if(p.size() == 0)
		throw std::runtime_error("Attempt to evaluate empty polynomial.");
	T value = 0;
	for(size_t i = 0; i < p.size(); ++i)
	{
		value = try_fma(argpow, p[i], value);
		argpow *= arg;
	}
	return value;
}

	
template <class T>
std::pair<T, T> root_bounds(const poly_t<T> & p)
{
	// from samuelusons inequality: 
	// https://en.wikipedia.org/wiki/Properties_of_polynomial_roots#Bounds_on_.28complex.29_polynomial_roots
	size_t n = p.degree();
	T midpoint = -p[n - 1] / ((big_float_t)n * p[n]);
	T span = ((n - 1) * boost::multiprecision::sqrt((p[n - 1] * p[n - 1]) - ((2 * n) * p[n] * p[n - 2])/ (n - 1))) / (n * p[n]);
	// 'lower bound' and 'upper bound' on the possible location of the roots
	T lb = midpoint - span;
	T ub = midpoint + span;
	if(lb > ub)
		std::swap(lb, ub);
	return {lb, ub};
}
 

template <class Real, class T>
poly_t<Real> as_real_poly(const poly_t<T> & p)
{
	poly_t<Real> rp;
	rp.data().resize(p.size());
	std::transform(p.data().begin(), p.data().end(), rp.data().begin(), [](const auto & arg){ return boost::rational_cast<Real>(arg);});
	return rp;
}

#endif /* POLYNOMIALS_H */
