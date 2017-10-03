#include "big_numbers.h"
#include <boost/functional/hash.hpp>
#include <vector>
#include <unordered_map>


using vec2_t = std::pair<size_t, size_t>;
struct DualHash
{
	boost::hash<size_t> hasher;
	size_t operator()(const vec2_t & pair) const
	{
		size_t seed = hasher(pair.first);
		boost::hash_combine(seed, pair.second);
		return seed;
	}
};
// arbitrary recursion limit on binomial coefficient calculation
static constexpr const size_t binom_coeff_max_recurse = 5;
using binom_cache_t = std::unordered_map<vec2_t, big_integer_t, DualHash>;


big_integer_t factorial(size_t n)
{
	std::vector<big_integer_t> cache{1, 1, 2, 6};
	if(cache.size() > n)
		return cache[n];
	cache.reserve(n);
	while(cache.size() <= n)
		cache.push_back(cache.back() * cache.size());
	return cache.back();
}

static big_integer_t binom_coeff_recurse(size_t n, size_t k, binom_cache_t & cache, size_t nrecurse)
{
	if((not k) or (n == k))
		return 1;
	auto & cached = cache[{n, k}];
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

big_integer_t binom_coeff(size_t n, size_t k)
{
	static binom_cache_t cache;
	return binom_coeff_recurse(n, k, cache, 0);
}
