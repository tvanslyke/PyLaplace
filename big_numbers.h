#ifndef BIG_NUMBERS_H
#define BIG_NUMBERS_H
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/rational.hpp>
#include <type_traits>
namespace{ namespace mp = boost::multiprecision; }
using big_integer_t = boost::multiprecision::cpp_int;
using big_float_t = mp::number<mp::backends::cpp_bin_float<256, mp::backends::digit_base_2, void, boost::int16_t, -1022, 1023>, mp::et_off>;
using rat_t = boost::rational<big_integer_t>;
big_integer_t factorial(size_t n);
big_integer_t binom_coeff(size_t n, size_t k);
template <class T, class U, class V>
auto try_fma(T a, U b, V c)
{
	constexpr bool a_good = std::is_floating_point_v<T> or std::is_integral_v<T>;
	constexpr bool b_good = std::is_floating_point_v<U> or std::is_integral_v<U>;
	constexpr bool c_good = std::is_floating_point_v<V> or std::is_integral_v<V>;
	if constexpr(a_good and b_good and c_good)
	{
		return std::fma(a, b, c);
	}
	else
	{
		return a * b + c;
	}
}


#endif /* BIG_NUMBERS_H */
