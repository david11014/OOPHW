#pragma once

#ifndef __myFloatingPoint_H__
#define __myFloatingPoint_H__

//////////////////////////////////////////////////////////////////////////
template <size_t size>
class TypeWithSize
{
public:
	typedef void Int;
};

template <>
class TypeWithSize<4>
{
public:
	typedef int Int;
	typedef unsigned int UInt;
};

template <>
class TypeWithSize<8> {
public:
#if WIN32
	typedef __int64 Int;
	typedef unsigned __int64 UInt;
#else
	typedef long long Int;  // NOLINT
	typedef unsigned long long UInt;  // NOLINT
#endif
};

template <typename RawType, size_t Ulps>
class FloatingPoint
{
public:
	// Defines the unsigned integer type that has the same size as the
	// floating point number.
	typedef typename TypeWithSize<sizeof(RawType)>::UInt Bits;
	typedef FloatingPoint<RawType, Ulps> Type;
	// Constants.

	// # of bits in a number.
	static const size_t kBitCount = 8 * sizeof(RawType);

	// # of fraction bits in a number.
	static const size_t kFractionBitCount =
		std::numeric_limits<RawType>::digits - 1;

	// # of exponent bits in a number.
	static const size_t kExponentBitCount = kBitCount - 1 - kFractionBitCount;

	// The mask for the sign bit.
	static const Bits kSignBitMask = static_cast<Bits>(1) << (kBitCount - 1);

	// The mask for the fraction bits.
	static const Bits kFractionBitMask =
		~static_cast<Bits>(0) >> (kExponentBitCount + 1);

	// The mask for the exponent bits.
	static const Bits kExponentBitMask = ~(kSignBitMask | kFractionBitMask);

	// How many ULP's (Units in the Last Place) we want to tolerate when
	// comparing two numbers.  The larger the value, the more error we
	// allow.  A 0 value means that two numbers must be exactly the same
	// to be considered equal.
	//
	// The maximum error of a single floating-point operation is 0.5
	// units in the last place.  On Intel CPU's, all floating-point
	// calculations are done with 80-bit precision, while double has 64
	// bits.  Therefore, 4 should be enough for ordinary use.
	//
	// See the following article for more details on ULP:
	// http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm.
	static const size_t kMaxUlps = Ulps;

	// Constructs a FloatingPoint from a raw floating-point number.
	//
	// On an Intel CPU, passing a non-normalized NAN (Not a Number)
	// around may change its bits, although the new value is guaranteed
	// to be also a NAN.  Therefore, don't expect this constructor to
	// preserve the bits in x when x is a NAN.
	FloatingPoint(const RawType& x)
	{
		u_.value_ = x;
	}

	// Static methods

	// Reinterprets a bit pattern as a floating-point number.
	//
	// This function is needed to test the AlmostEquals() method.
	static RawType ReinterpretBits(const Bits bits)
	{
		FloatingPoint fp(0);
		fp.u_.bits_ = bits;
		return fp.u_.value_;
	}

	// Returns the floating-point number that represent positive infinity.
	static RawType Infinity()
	{
		return ReinterpretBits(kExponentBitMask);
	}

	// Non-static methods

	// Returns the bits that represents this number.
	const Bits &bits() const
	{
		return u_.bits_;
	}

	// Returns the exponent bits of this number.
	Bits exponent_bits() const
	{
		return kExponentBitMask & u_.bits_;
	}

	// Returns the fraction bits of this number.
	Bits fraction_bits() const
	{
		return kFractionBitMask & u_.bits_;
	}

	// Returns the sign bit of this number.
	Bits sign_bit() const
	{
		return kSignBitMask & u_.bits_;
	}

	// Returns true iff this is NAN (not a number).
	bool is_nan() const
	{
		// It's a NAN if the exponent bits are all ones and the fraction
		// bits are not entirely zeros.
		return (exponent_bits() == kExponentBitMask) && (fraction_bits() != 0);
	}

	// Returns true iff this number is at most kMaxUlps ULP's away from
	// rhs.  In particular, this function:
	//
	//   - returns false if either number is (or both are) NAN.
	//   - treats really large numbers as almost equal to infinity.
	//   - thinks +0.0 and -0.0 are 0 DLP's apart.
	bool AlmostEquals(const FloatingPoint& rhs) const
	{
		// The IEEE standard says that any comparison operation involving
		// a NAN must return false.
		if (is_nan() || rhs.is_nan())
		{
			return false;
		}

		return DistanceBetweenSignAndMagnitudeNumbers(u_.bits_, rhs.u_.bits_) <= kMaxUlps;
	}

	bool operator==(const FloatingPoint& rhs) const
	{
		return u_.value_ == rhs.u_.value_ || AlmostEquals(rhs);
	}

	bool operator==(const RawType& rhs) const
	{
		return AlmostEquals(FloatingPoint(rhs));
	}

	bool operator<(const FloatingPoint& rhs) const
	{
		return u_.value_ < rhs.u_.value_ && !AlmostEquals(rhs);
	}

	bool operator>(const FloatingPoint& rhs) const
	{
		return u_.value_ > rhs.u_.value_ && !AlmostEquals(rhs);
	}

	bool operator<=(const FloatingPoint& rhs) const
	{
		return u_.value_ <= rhs.u_.value_ || AlmostEquals(rhs);
	}

	bool operator>=(const FloatingPoint& rhs) const
	{
		return u_.value_ >= rhs.u_.value_ || AlmostEquals(rhs);
	}
	Type operator*(const Type& rhs) const
	{
		return u_.value_*rhs.u_.value_;
	}
	Type operator/(const Type& rhs) const
	{
		return u_.value_/rhs.u_.value_;
	}
	operator double() const
	{
		return u_.value_;
	}
private:
	// The data type used to store the actual floating-point number.
	union FloatingPointUnion
	{
		RawType value_;  // The raw floating-point number.
		Bits bits_;      // The bits that represent the number.
	};

	// Converts an integer from the sign-and-magnitude representation to
	// the biased representation.  More precisely, let N be 2 to the
	// power of (kBitCount - 1), an integer x is represented by the
	// unsigned number x + N.
	//
	// For instance,
	//
	//   -N + 1 (the most negative number representable using
	//          sign-and-magnitude) is represented by 1;
	//   0      is represented by N; and
	//   N - 1  (the biggest number representable using
	//          sign-and-magnitude) is represented by 2N - 1.
	//
	// Read http://en.wikipedia.org/wiki/Signed_number_representations
	// for more details on signed number representations.
	static Bits SignAndMagnitudeToBiased(const Bits &sam)
	{
		if (kSignBitMask & sam) {
			// sam represents a negative number.
			return ~sam + 1;
		}
		else {
			// sam represents a positive number.
			return kSignBitMask | sam;
		}
	}

	// Given two numbers in the sign-and-magnitude representation,
	// returns the distance between them as an unsigned number.
	static Bits DistanceBetweenSignAndMagnitudeNumbers(const Bits &sam1, const Bits &sam2)
	{
		const Bits biased1 = SignAndMagnitudeToBiased(sam1);
		const Bits biased2 = SignAndMagnitudeToBiased(sam2);
		return (biased1 >= biased2) ? (biased1 - biased2) : (biased2 - biased1);
	}

	FloatingPointUnion u_;
};

typedef FloatingPoint<float, 6> myFloat;
typedef FloatingPoint<double, 6> myDouble;

template<size_t N>
using myFloatt = FloatingPoint<float, N>;

template<size_t N>
using myDoublet = FloatingPoint<double, N>;
//////////////////////////////////////////////////////////////////////////

#endif