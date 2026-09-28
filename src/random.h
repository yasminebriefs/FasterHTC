#pragma once

#include <optional>
#include <random>

namespace fasterhtc {

class RandomTool {
public:
	RandomTool (uint64_t, std::optional<int64_t>, bool (*)(int64_t));

	int64_t getRandom (int64_t, int64_t);

	int64_t getRandomSmallerPrime ();

	int64_t getPrime () const noexcept;

private:
	std::mt19937_64 rng;
	const int64_t prime;
};

}
