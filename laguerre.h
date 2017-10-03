#ifndef LAGUERRE_H
#define LAGUERRE_H
#include <algorithm>
#include "polynomials.h"
#include "precomputed_laguerre_roots.h"


template <class DestIt>
big_integer_t laguerre_poly(DestIt begin, DestIt end)
{
	// sum{  n! / ((k!)^2 * (n - k)!) 
	size_t n = std::distance(begin, end);
	int_least8_t sign = 1;
	big_integer_t k_factorial = 1;
	big_integer_t n_factorial = factorial(n);
	big_integer_t n_minus_k_factorial;
	big_integer_t coeff_term = n_factorial;
	for(size_t k = 0; k < n;)
	{
		// TODO: optimize this with binom_coeff(n, k) and factorial(k)
		*begin++ = coeff_term * sign;
		sign *= -1;
		++k;
		coeff_term *= n - k;	
		coeff_term /= (k * k);
	}
	assert(std::distance(begin, end) == 0);
	return n_factorial;
}

void _compute_laguerre_roots(size_t n, std::vector<std::vector<big_float_t>>& cache);

template <class DestIt>
DestIt laguerre_roots(size_t n, DestIt dest)
{
	static std::vector<std::vector<big_float_t>>& cache = precomputed_laguerre_roots();
	if(n >= cache.size())
		_compute_laguerre_roots(n, cache);
	std::transform(cache[n].begin(), cache[n].end(), dest, [](const auto & v){ return (long double)(v);});
	return dest;
}
inline const std::vector<big_float_t> & cached_laguerre_roots(size_t n)
{
	static const std::vector<std::vector<big_float_t>>& cache = precomputed_laguerre_roots();
	return cache[n];
}


template <class WeightType, class PolyFunc, class It, class DestIt>
DestIt laguerre_weights(PolyFunc pf, It roots_begin, It roots_end, DestIt dest)
{
	big_integer_t n_plus_one_sqrd = std::distance(roots_begin, roots_end);
	++n_plus_one_sqrd;
	n_plus_one_sqrd *= n_plus_one_sqrd;
	WeightType denom_coeff = static_cast<double>(n_plus_one_sqrd);
	WeightType eval_result = 0;
	auto transform_func = [&](const auto & root)
	{ 
		eval_result = pf(root);
		return static_cast<WeightType>(root / (denom_coeff * (eval_result * eval_result)));
	};
	return std::transform(roots_begin, roots_end, dest, transform_func);
}


#endif /* LAGUERRE_H */
