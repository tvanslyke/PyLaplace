#include <cmath>
#include <vector>
#include <unordered_map>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/functional/hash.hpp>
#include <boost/rational.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <array>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <complex>
namespace mp = boost::multiprecision;
namespace math = boost::math::tools;

mp::cpp_int factorial(size_t n)
{
	std::vector<mp::cpp_int> cache{1, 1, 2, 6};
	if(cache.size() > n)
		return cache[n];
	else
	{
		while(not (cache.size() > n))
			cache.push_back(cache.back() * cache.size());
		return cache.back();
	}
}
namespace {
	using vec2_t = std::array<size_t, 2>;
	struct DualHash
	{
		boost::hash<size_t> hasher;
		size_t operator()(const std::array<size_t, 2> & pair) const
		{
			size_t seed = hasher(pair[0]);
			boost::hash_combine(seed, hasher(pair[1]));
			return seed;
		}
	};
	// arbitrary recursion limit on binomial coefficient calculation
	static constexpr const size_t binom_coeff_max_recurse = 5;
	using binom_cache_t = std::unordered_map<vec2_t, mp::cpp_int, DualHash>;
	using rat_t = boost::rational<mp::cpp_int>;
	using poly_t = math::polynomial<rat_t>;
	mp::cpp_int binom_coeff_recurse(size_t n, size_t k, binom_cache_t & cache, size_t nrecurse)
	{
		if((not k) or (n == k))
			return 1;
		auto & cached = cache[vec2_t{n, k}];
		if(not cached)
		{
			if(nrecurse < binom_coeff_max_recurse)
				cached = binom_coeff_recurse(n - 1, k - 1, cache, nrecurse + 1) + 
			       		 binom_coeff_recurse(n - 1, k, cache, nrecurse + 1);
			else
				cached = factorial(n) / (factorial(k) * factorial(n - k));
		}
		return cached;
	}
} /* anonymous namespace */

struct PolySolverFunc
{
	using triple_t = boost::math::tuple<double, double, double>;
	template <class ... Args>
	triple_t operator()(double x) const
	{
		const auto & p = poly_;
		if(p.size() == 1)
			return triple_t{p[0], 0.0, 0.0};
		else if(p.size() == 2)
			return triple_t{p[0] + p[1] * x, p[1], 0.0};
		double ddf = p[2] / 2;
		double df = std::fma((x / 2), p[2], p[1]);
		double f =  std::fma(p[1], x, std::fma(p[2], x * x, p[0])); // p[0] + x * p[1] + x * x * p[2]
		size_t n = 3;
		double x1 = x;
		double x2 = x * x;
		double x3 = x * x * x;
		while(n < p.size())
		{
			ddf = std::fma((x1 / n) / (n - 1), p[n], ddf);
			df = std::fma(x2 / n, p[n], df);
			f = std::fma(x3, p[n], f);
			x1 = x2;
			x2 = x3;
			x3 *= x;
			++n;
		}
		return triple_t{f, df, ddf};
	}

	math::polynomial<double> poly_;
};

template <class It, class T>
auto evaluate_polynomial(It coeffs_begin, It coeffs_end, T arg)
{
	// TODO: bounds checking on distance(coeffs_begin, coeffs_end);
	auto tmp = *coeffs_begin++;
	if(not (coeffs_begin != coeffs_end))
	{
		return tmp;	
	}
	else
	{
		using value_t = std::conditional_t<std::is_integral_v<T>, mp::cpp_int, T>;
		value_t argpow = arg;
		auto result = tmp + argpow * (*coeffs_begin++);
		while(coeffs_begin != coeffs_end)
			result += (*coeffs_begin++) * (argpow *= arg);
		return result;
	}
}
template <class T>
auto evaluate_polynomial(const math::polynomial<T> & p, T arg)
{
	T argpow = 1;
	T value = 0;
	for(size_t i = 0; i < p.size(); ++i)
	{
		if constexpr(std::is_floating_point_v<T>){
			value = std::fma(argpow, p[i], value);	
		}else{
			value += argpow * p[i];
		}
		argpow *= arg;
	}
	return value;
}

template <class Real = double>
math::polynomial<Real> as_real_poly(const poly_t & p)
{
	math::polynomial<Real> rp;
	rp.data().resize(p.size());
	std::transform(p.data().begin(), p.data().end(), rp.data().begin(), [](const auto & arg){ return boost::rational_cast<Real>(arg);});
	return rp;
}

math::polynomial<rat_t> poly_diff(const math::polynomial<rat_t> & poly)
{
	math::polynomial<rat_t> diff(poly.data().begin() + 1, poly.data().end());
	for(size_t i = 0; i < diff.size(); ++i)
		diff[i] /= (i + 1);
	return diff;	
}

math::polynomial<double> poly_diff(const math::polynomial<double> & poly)
{
	math::polynomial<double> diff(poly.data().begin() + 1, poly.data().end());
	for(size_t i = 0; i < diff.size(); ++i)
		diff[i] /= (i + 1);
	return diff;	
}

auto poly_solver_func(const math::polynomial<double>& p)
{
	auto psf = PolySolverFunc();
	psf.poly_ = p;
	return psf;
}
auto poly_solver_func(poly_t rat_poly)
{
	return poly_solver_func(as_real_poly(rat_poly));
}


mp::cpp_int binom_coeff(size_t n, size_t k)
{
	static binom_cache_t cache;
	return binom_coeff_recurse(n, k, cache, 0);
}

const poly_t & recursive_laguerre_polynomial(size_t order)
{
	static const poly_t one_minus_x({1, -1});
	static std::vector<poly_t> cache{{1}, {1, -1}};
	if(cache.size() > order)
		return cache[order];
	const rat_t k(order);
	poly_t value = ((2 * (k - 1) + one_minus_x) * recursive_laguerre_polynomial(order - 1) - 
		       (k - 1) * recursive_laguerre_polynomial(order - 2)) / k;
	cache.push_back(value);
	assert(cache.size() == order + 1);
	return cache.back();
}


template <class DestIt>
mp::cpp_int laguerre_poly(DestIt begin, DestIt end)
{
	// sum{  n! / ((k!)^2 * (n - k)!) 
	size_t n = std::distance(begin, end);
	int_least8_t sign = 1;
	mp::cpp_int k_factorial = 1;
	mp::cpp_int n_factorial = factorial(n);
	mp::cpp_int n_minus_k_factorial;
	mp::cpp_int coeff_term = n_factorial;
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



std::pair<double, double> root_bounds(const math::polynomial<double>& p)
{
	// from samuelusons inequality: 
	// https://en.wikipedia.org/wiki/Properties_of_polynomial_roots#Bounds_on_.28complex.29_polynomial_roots
	size_t n = p.degree();
	double midpoint = -1 * p[n - 1] / (n * p[n]);
	double span = ((n - 1) * std::sqrt((p[n - 1] * p[n - 1]) - ((2 * n) * p[n] * p[n - 2])/ (n - 1))) / (n * p[n]);
	// 'lower bound' and 'upper bound' on the possible location of the roots
	double lb = midpoint - span;
	double ub = midpoint + span;
	if(lb > ub)
		std::swap(lb, ub);
	return {lb, ub};
}


void compute_laguerre_roots(size_t n, std::vector<std::vector<double>>& cache)
{
	if(cache.size() < n)
		compute_laguerre_roots(n - 1, cache);
	else if(cache.size() > n)
		return;
	
	poly_t rat_poly;
	rat_poly.data().resize(n + 1);
	laguerre_poly(rat_poly.data().begin(), rat_poly.data().end());
	auto p = as_real_poly<double>(rat_poly);
	
	auto obj_func = [&](double v)
	{
		return evaluate_polynomial(p, v);
	};
	auto tol_func = [&](double left, double right)
	{
		int cls = std::fpclassify(right - left);
		if((cls == FP_ZERO) or (cls == FP_SUBNORMAL))
			return true;
		cls = std::fpclassify(obj_func(left + (right - left) / 2));
		return ((cls == FP_ZERO) or (cls == FP_SUBNORMAL));
	};
	auto pair_mean = [](const std::pair<double, double> & pair)
	{
		return std::get<0>(pair) + (std::get<1>(pair) - std::get<0>(pair)) / 2;
	};
	rat_poly.data().resize(0);
	rat_poly.data().shrink_to_fit();
	auto [lb, ub] = root_bounds(p);
	lb = std::max(lb, 0.0);
	cache.push_back(std::vector<double>(n));
	auto prev_roots = cache[cache.size() - 2];
	auto dest = cache.back().begin();
	auto prev_iter = prev_roots.begin();
	*dest++ = pair_mean(math::bisect(obj_func, lb, *prev_iter++, tol_func));
	while(prev_iter < prev_roots.end())
	{
		*dest++ = pair_mean(math::bisect(obj_func, *(prev_iter - 1), *prev_iter, tol_func));
		++prev_iter;
	}
	*dest++ = pair_mean(math::bisect(obj_func, prev_roots.back(), ub, tol_func));
}
template <class DestIt>
DestIt laguerre_roots(size_t n, DestIt dest)
{
	static std::vector<std::vector<double>> cache(
		[](){
			std::vector<std::vector<double>> v;
			v.push_back(std::vector<double>());
			v.push_back(std::vector<double>{1.0});
			v.push_back(std::vector<double>{0.585786437627, 3.414213562373});
			return v;
		}()
	);
	if(n >= cache.size())
		compute_laguerre_roots(n, cache);
	std::copy(cache[n].begin(), cache[n].end(), dest);
	return dest;
}
template <class DestIt>
DestIt laguerre_roots(const poly_t & p, DestIt dest)
{
	return laguerre_roots(as_real_poly(p), dest);
}
template <class InputIt, class DestIt>
DestIt laguerre_roots(InputIt begin, InputIt end, DestIt dest)
{
	math::polynomial<double> poly;
	poly.data().assign(begin, end);
	return laguerre_roots(poly, dest);
}


auto poly_pow(const poly_t & p, size_t i)
{
	if(i == 0)
		return poly_t{1};
	else if(i == 1)
		return p;
	else
	{
		auto result = poly_pow(p, i / 2);
		result *= result;
		if(i % 2 == 0)
			result *= p;
		return result;
	}

}



template <class T, class U>
auto evaluate_polynomial(const math::polynomial<T> & p, U && arg)
{
	return evaluate_polynomial(p.data().begin(), p.data().end(), std::forward<U>(arg));
}

template <class PolyFunc, class It, class DestIt>
DestIt laguerre_weights(PolyFunc pf, It roots_begin, It roots_end, DestIt dest)
{
	mp::cpp_int n_plus_one_sqrd = std::distance(roots_begin, roots_end);
	++n_plus_one_sqrd;
	n_plus_one_sqrd *= n_plus_one_sqrd;
	double denom_coeff = static_cast<double>(n_plus_one_sqrd);
	double eval_result = 0;
	auto transform_func = [&](const auto & root)
	{ 
		eval_result = pf(root);
		return root / (denom_coeff * (eval_result * eval_result));
	};
	return std::transform(roots_begin, roots_end, dest, transform_func);
}



template <class F, class WeightIt, class RootIt>
auto evaluate_gauss_laguerre_quadrature(F func, WeightIt wbegin, WeightIt wend, RootIt rbegin, RootIt rend)
{
	auto wdist = std::distance(wbegin, wend);
	if(wdist < 1)
		throw std::invalid_argument("Must have at least one weight and one root in "
					    "'evaluate_gauss_laguerre_quadrature()'.");
	if(wdist != std::distance(rbegin, rend))
		throw std::invalid_argument("Must have same number of weights and roots in "
	        		            "'evaluate_gauss_laguerre_quadrature()'.");

	auto result = *wbegin * func(*rbegin);
	++wbegin;
	++rbegin;
	while(wbegin != wend)
	{
		if constexpr(std::is_floating_point_v<decltype(result)>)
		{
			result = std::fma(*wbegin, func(*rbegin), result);
		}
		else
		{
			result += (*wbegin) * func(*rbegin);
		}
		++wbegin;
		++rbegin;
	}
	return result;
}

	

template <class It>
auto poly_func(It begin, It end)
{
	using value_t = typename std::iterator_traits<It>::value_type;
	const std::vector<value_t> coeffs;
	return [=](const auto & arg)
	{
		return evaluate_polynomial(coeffs.begin(), coeffs.end(), arg);
	};
}



template <class F>
auto laplace_transform(F func, size_t order)
{
	// g(t, s) = f(t) * exp(-1 * s * t) * dt
	// g(t, s) = (f(t) * exp(-1 * (s - 1) * t)) * exp(-t) * dt
	// F = integral(f(t) * exp(-t), t)
	// integral(g(t, s), t) = F
	std::vector<double> roots(order);
	std::vector<double> feval_roots(order);
	std::vector<double> weights(order);
	poly_t tmp;
	tmp.data().resize(order + 2);
	laguerre_roots(order, roots.begin());
	auto div = laguerre_poly(tmp.data().begin(), tmp.data().end());
	tmp /= div;
	auto poly = as_real_poly(tmp);
	tmp.data().resize(0);
	tmp.data().shrink_to_fit();
	auto poly_func = [=](double d)
	{
		return evaluate_polynomial(poly, d);
	};

	laguerre_weights(poly_func, roots.begin(), roots.end(), weights.begin());
	std::transform(roots.begin(), roots.end(), feval_roots.begin(), func);
	for(size_t i = 0; i < weights.size(); ++i)
	{
		weights[i] *= func(roots[i]);
	}
	return [=](std::complex<double>  s)
	{
		std::complex<double> result{0.0, 0.0};
		s = std::complex<double>{1, 0} - s;
		for(size_t i = 0; i < roots.size(); ++i)
			result += weights[i] * std::exp(roots[i] * s);
		return result;
	};
}


#include <iomanip>

int main()
{
	size_t n = 20;
	auto func = [](double t){ return std::exp(-2 * t); };
	auto tfunc = laplace_transform(func, n);
	auto exact = [](double s){ return 1 / (s + 2); };
	for(double s = 0.2; s < 30.0; s += 0.1)
	{
		std::cout << tfunc(s).real() << " " << exact(s) << std::endl;
	}


	return 0;
}


