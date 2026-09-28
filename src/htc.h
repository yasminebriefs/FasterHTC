#pragma once

#include "algebra.h"
#include "random.h"

#include <optional>
#include <variant>
#include <vector>

namespace fasterhtc {

struct IdentificationResult {
	std::vector<std::pair<size_t, size_t>> identification;
	uint64_t seed;
	int64_t prime;

	IdentificationResult (uint64_t seed_, int64_t prime_) : identification(0), seed(seed_), prime(prime_) {}

	void addIdentifiedEdge (size_t i, size_t j) {
		this->identification.emplace_back(i, j);
	}

	void addIdentifiedEdgeTranslated (size_t i, size_t j, const std::vector<size_t> &translation) {
		this->identification.emplace_back(translation[i], translation[j]);
	}
};

IdentificationResult htc (size_t n, const std::vector<std::pair<size_t, size_t>> &bidirected, const std::vector<std::pair<size_t, size_t>> &directed, std::optional<uint64_t> seed, std::optional<int64_t> prime, bool (*)(int64_t));

}
