#pragma once
#include <vector>
#include <string>
#include "bigint/limb.h"

namespace bi
{
	struct BigInt
	{
		// little-endian limbs (limbs[0] is least significant)
		std::vector<limb_t> limbs;
		BigInt() = default;
		explicit BigInt(std::size_t n, limb_t v = 0) : limbs(n, v) {}
		explicit BigInt(limb_t v)
		{
			if (v)
				limbs.push_back(v);
		}

		// static BigInt from_u64(std::uint64_t x) { return BigInt(x); }
		static BigInt from_u64(std::uint64_t x)
		{
			BigInt r;
			if (x)
				r.limbs.push_back((limb_t)x);
			return r;
		}
		void normalize();
		std::string to_hex() const;
	};

	// basic operations
	void add_inplace(BigInt &a, const BigInt &b);
	BigInt add(const BigInt &a, const BigInt &b);
	BigInt mul_schoolbook(const BigInt &a, const BigInt &b);
	int cmp(const BigInt &a, const BigInt &b);
	bool sub_inplace(BigInt &a, const BigInt &b);
	std::uint64_t to_u64(const BigInt& x);
}
