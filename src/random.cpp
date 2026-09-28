#include "random.h"

#include <stdexcept>

namespace fasterhtc {

uint64_t deriveGeneralSeed (uint64_t seed) {
	std::mt19937_64 gen(seed);
	return gen();
}

uint64_t derivePrimeSeed (uint64_t seed) {
	std::mt19937_64 gen(seed);
	gen();
	return gen();
}

int64_t generateCandidate (std::mt19937_64 &gen) {
	return std::uniform_int_distribution<int64_t>(1ll << 58, (1ll << 59) - 1)(gen) | 1ll;
}

int64_t generatePrime (uint64_t prime_seed, bool (*isProbablyPrime)(int64_t)) {
	std::mt19937_64 gen(prime_seed);

	int64_t candidate = generateCandidate(gen);
	while (!isProbablyPrime(candidate)) {
		candidate = generateCandidate(gen);
	}
	return candidate;
}

RandomTool::RandomTool (uint64_t seed, std::optional<int64_t> prime_, bool (*isProbablyPrime)(int64_t)) : rng(deriveGeneralSeed(seed)), prime(prime_.has_value() ? prime_.value() : generatePrime(derivePrimeSeed(seed), isProbablyPrime)) {}

int64_t RandomTool::getRandom (int64_t l, int64_t r) {
	if (l > r) {
		throw std::logic_error("getRandom: Call with invalid bounds. Please report this error.");
	}
	return std::uniform_int_distribution<int64_t>(l, r)(this->rng);
}

int64_t RandomTool::getRandomSmallerPrime () {
	return this->getRandom(1, prime - 1);
}

int64_t RandomTool::getPrime () const noexcept {
	return this->prime;
}

}
