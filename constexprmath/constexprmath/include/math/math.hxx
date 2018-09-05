#pragma once

#include <type_traits>
#include <limits>
#include <array>

#define LN2 0.69314718055994530941723212145818
#define PI	3.1415926535897932384626433832795

using integral = std::int64_t;

namespace math
{
	template <typename T>
	static constexpr std::enable_if_t<std::is_arithmetic_v<T>, integral> log2i(T x) noexcept
	{
		integral count = 0, val = 1;

		while (T(val) <= x)
		{
			count++;
			val <<= 1;
		}

		return count - 1;
	}

	template <typename T>
	constexpr T max(const T& fst, const T& snd)
	{
		return fst > snd ? fst : snd;
	}

	template <typename T>
	constexpr T min(const T& fst, const T& snd)
	{
		return fst < snd ? fst : snd;
	}

	template <typename T>
	constexpr auto log(T x) noexcept
	{
		using RCV = std::remove_cv_t<std::remove_reference_t<T>>;

		if constexpr (std::is_integral_v<RCV>)
		{
			return log(static_cast<T>(x));
		}
		else
		{
			if (x < RCV(0))
			{
				return std::numeric_limits<RCV>::quiet_NaN();
			}

			if (x == RCV(0))
			{
				return -std::numeric_limits<RCV>::infinity();
			}

			if (x == RCV(1))
			{
				return T(0);
			}

			if (x == std::numeric_limits<RCV>::infinity())
			{
				return std::numeric_limits<RCV>::infinity();
			}

			if (x < RCV(1))
			{
				return -log(RCV(1) / x);
			}

			integral l2 = log2i(x);

			const auto new_x = x / RCV(1 << l2);


			const auto zratio = (new_x - 1) / (new_x + 1);
			const auto zratio2 = zratio * zratio;
			const auto zratio3 = zratio2 * zratio;
			const auto zratio4 = zratio2 * zratio2;
			const auto zratio5 = zratio4 * zratio;
			const auto zratio6 = zratio4 * zratio2;
			const auto zratio7 = zratio5 * zratio2;
			const auto zratio8 = zratio4 * zratio4;
			const auto zratio9 = zratio5 * zratio4;

			return 2 * (
				zratio
				+ zratio3 / 3
				+ zratio5 / 5
				+ zratio7 / 7
				+ zratio5 * zratio4 / 9
				+ zratio7 * zratio4 / 11
				+ zratio7 * zratio6 / 13
				+ zratio7 * zratio8 / 15
				+ zratio9 * zratio8 / 17
				+ zratio9 * zratio8 * zratio2 / 19
				+ zratio9 * zratio8 * zratio4 / 21
				+ zratio9 * zratio8 * zratio6 / 23
				) + l2 * T(LN2);
		}
	}

	template <typename T>
	constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> abs(T x) noexcept
	{
		return x < 0 ? -x : x;
	}

	template <typename T>
	static constexpr std::enable_if_t<std::is_floating_point_v<T>, bool> satisfies_min_error(T x, T val) noexcept
	{
		if constexpr (std::is_same_v<T, float>)
		{
			return abs(x - val) < 0.000000000000001F;
		}
		if constexpr (std::is_same_v<T, double>)
		{
			return abs(x - val) < 0.000000000000001;
		}
		if constexpr (std::is_same_v<T, long double>)
		{
			return abs(x - val) < 0.000000000000005L;
		}
	}

	template <typename T>
	constexpr auto sqrt(T x) noexcept
	{
		if constexpr (std::is_integral_v<T>)
		{
			return sqrt(static_cast<double>(x));
		}
		else
		{
			if (x < T(0))
			{
				return std::numeric_limits<T>::quiet_NaN();
			}

			if (x == std::numeric_limits<T>::infinity())
			{
				return std::numeric_limits<T>::infinity();
			}

			T result = T(1);

			for (int i = 0; i < 20; ++i)
			{
				const auto r2 = result * result;
				result -= (2 * r2 * result - 2 * x * result) / (3 * r2 + x);
			}

			return result;
		}
	}

	/**
	  *	O(1) Exponential function that can function both as a runtime and compile-time means
	  * of estimating a highly accurate value of e^x, although this accuracy becomes reduced
	  * at exceedingly higher values of x.
	  */
	template <typename T>
	constexpr auto exp(T x) noexcept
	{

		if (x == std::numeric_limits<T>::infinity())
		{
			return std::numeric_limits<T>::infinity();
		}

		if (x == -std::numeric_limits<T>::infinity())
		{
			return T(0);
		}

		if (x < T(0))
		{
			return 1 / exp(-x);
		}

		if (x == T(0))
		{
			return T(1);
		}

		if (x > T(0) && x < T(1))
		{
			const auto x1 = T(1 / 2.0) * x;
			const auto x2 = T(1 / 4.0) * x * x1;
			const auto x3 = T(1 / 6.0) * x * x2;
			const auto x4 = T(1 / 8.0) * x * x3;
			const auto x5 = T(1 / 10.0) * x * x4;
			const auto x6 = T(1 / 12.0) * x * x5;
			const auto x7 = T(1 / 14.0) * x * x6;
			const auto x8 = T(1 / 16.0) * x * x7;
			const auto x9 = T(1 / 18.0) * x * x8;
			const auto x10 = T(1 / 20.0) * x * x9;
			const auto x11 = T(1 / 22.0) * x * x10;

			return (1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11)
				/ (1 - x1 + x2 - x3 + x4 - x5 + x6 - x7 + x8 - x9 + x10 - x11);
		}

		else
		{
			constexpr integral sz = sizeof(T) * 8 - 2;
			integral k = static_cast<integral> (x / T(LN2));
			const T rem = x - T(k) * T(LN2);
			T rexp = exp(rem);

			while (k > 0)
			{
				rexp *= integral(1) << (min(k, sz));
				k -= sz;
			}

			return rexp;
		}
		
	}

	template <typename T>
	constexpr auto sigmoid(const T &x) noexcept
	{
		return T(1) / (T(1) + exp(-x));
	}

	constexpr double E = exp(1.0);

	template <typename T>
	constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> pow(T x, T y) noexcept
	{
		if (x == T(0))
		{
			return y == T(0) ? 1 : 0;
		}

		return x == E ? exp(y) : exp(y * log(x));
	}

	template <typename T>
	constexpr T round(const T &x) noexcept
	{
		return abs(static_cast<integral>(x) - x) > abs(static_cast<integral>(x) + 1 - x)
			? T(static_cast<integral>(x) + 1)
			: T(static_cast<integral>(x));
	}
}