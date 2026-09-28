#pragma once

#include "random.h"

#include <vector>

namespace fasterhtc {

int64_t multiply (int64_t, int64_t, int64_t) noexcept;

int64_t inverse (int64_t, int64_t) noexcept;

int64_t fastExp (int64_t, int64_t, int64_t) noexcept;

std::vector<std::vector<int64_t>> matrixMultiply (const std::vector<std::vector<int64_t>> &, const std::vector<std::vector<int64_t>> &, int64_t);

std::vector<std::vector<int64_t>> matrixInverse (std::vector<std::vector<int64_t>>, int64_t);

}
