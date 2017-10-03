#include <cmath>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <complex>
#include "laguerre.h"

template <class T>
class LaplaceTransform
{
public:
	using real_t = T;
	using complex_t = std::complex<real_t>;
	using paired_values_t = std::vector<std::pair<real_t, real_t>>;
	LaplaceTransform() = default;
	LaplaceTransform(const std::vector<real_t>& roots, const std::vector<real_t>& weights):
		root_weight_pairs_(zip(roots, weights))
	{

	}
	complex_t operator()(complex_t s) const
	{
		complex_t result{0.0, 0.0};
		complex_t exp_term;
		s.real(1 - s.real());
		for(const auto& [root, weight]: root_weight_pairs_)
		{
			// split up computation to allow the use of fused multiply-add instructions
			exp_term = std::exp(root * s);
			result.real(try_fma(weight, exp_term.real(), result.real()));
			result.imag(try_fma(weight, exp_term.imag(), result.imag()));
		}
		return result;
	}
	complex_t operator()(real_t s) const
	{
		real_t result = 0.0;
		s = 1 - s;
		for(const auto & [root, weight]: root_weight_pairs_)
			result = try_fma(weight, std::exp(root * s), result);
		return complex_t{result, 0.0};
	}
private:
	static paired_values_t zip(const std::vector<real_t>& left, const std::vector<real_t>& right)
	{
		if(left.size() != right.size())
			throw std::invalid_argument("Attempt to create LaplaceTransform object with different number of roots and weights.");
		paired_values_t zipped(left.size());
		for(size_t i = 0; i < left.size(); ++i)
		{
			zipped[i] = std::make_pair(left[i], right[i]);
		}
		return zipped;
	}
	const paired_values_t root_weight_pairs_;
};


template <class F>
auto laplace_transform(F func, size_t order)
{

	std::vector<double> roots(order);
	std::vector<double> feval_roots(order);
	std::vector<double> weights(order);
	poly_t<rat_t> tmp;
	tmp.data().resize(order + 2);
	laguerre_roots(order, roots.begin());
	auto div = laguerre_poly(tmp.data().begin(), tmp.data().end());
	tmp /= div;
	auto poly = as_real_poly<double>(tmp);
	tmp.data().resize(0);
	tmp.data().shrink_to_fit();
	auto poly_func = [=](double d)
	{
		return evaluate_polynomial(poly, d);
	};

	laguerre_weights<double>(poly_func, roots.cbegin(), roots.cend(), weights.begin());
	for(size_t i = 0; i < weights.size(); ++i)
	{
		weights[i] *= func(roots[i]);
	}
	return LaplaceTransform<double>(roots, weights);
}


