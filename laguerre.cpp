#include "laguerre.h"
#include <iostream>
namespace math = boost::math::tools;
void _compute_laguerre_roots(size_t n, std::vector<std::vector<big_float_t>>& cache)
{
	static const big_float_t tolerance_value = 0.0;
	if(cache.size() < n)
		_compute_laguerre_roots(n - 1, cache);
	else if(cache.size() > n)
		return;
	
	poly_t<rat_t> rational_poly;
	rational_poly.data().resize(n + 1);
	auto div = laguerre_poly(rational_poly.data().begin(), rational_poly.data().end());
	// rational_poly /= div;
	auto p = as_real_poly<big_float_t>(rational_poly);
	rational_poly.data().resize(0);
	rational_poly.data().shrink_to_fit();
	
	auto obj_func = [&](big_float_t v)
	{
		return evaluate_polynomial(p, v);
	};
	auto pair_mean = [](const std::pair<big_float_t, big_float_t> & pair)
	{
		return std::get<0>(pair) + (std::get<1>(pair) - std::get<0>(pair)) / 2;
	};
	auto tol_func = [&](big_float_t left, big_float_t right)
	{
		assert(right > left);
		return (right - left) < tolerance_value;
	};
	auto [lb, ub] = root_bounds(p);
	lb = std::max(lb, static_cast<big_float_t>(0.0));
	cache.emplace_back(n);
	const auto& prev_roots = cache[cache.size() - 2];
	auto dest = cache.back().begin();
	auto prev_iter = prev_roots.cbegin();
	*dest++ = pair_mean(math::bisect(obj_func, lb, *prev_iter++, tol_func));
	while(prev_iter < prev_roots.end())
	{
		*dest++ = pair_mean(math::bisect(obj_func, *(prev_iter - 1), *prev_iter, tol_func));
		++prev_iter;
	}
	*dest++ = pair_mean(math::bisect(obj_func, prev_roots.back(), ub, tol_func));
}
