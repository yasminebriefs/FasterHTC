#include "algebra.h"

#include <stdexcept>

namespace fasterhtc {

int64_t multiply (int64_t a, int64_t b, int64_t prime) noexcept {
#if defined(__SIZEOF_INT128__)

	__extension__ __int128 A = a, B = b;
	return (int64_t)(A * B % prime);

#else

	auto multiplyReduce = [](int64_t left, int64_t right, int64_t prime) {
		/* Assumes 0 <= left, right < 2^30 and 2^58 < prime < 2^59. Returns left * right % prime */
		int64_t ans = left * right;
		if (ans >= prime) ans -= prime;
		if (ans >= prime) ans -= prime;
		if (ans >= prime) ans -= prime;
		return ans;
	};

	auto double30Times = [](int64_t num, int64_t prime) {
		/* Assumes 0 <= num < prime. Returns num * 2^30 % prime */
		for (size_t i = 0; i < 30; i++) {
			num *= 2;
			if (num >= prime) num -= prime;
		}
		return num;
	};

	int64_t a0 = a & ((1 << 30) - 1), a1 = a >> 30;
	int64_t b0 = b & ((1 << 30) - 1), b1 = b >> 30;

	int64_t tmp = multiplyReduce(a0, b1, prime) + multiplyReduce(a1, b0, prime) + double30Times(multiplyReduce(a1, b1, prime), prime);
	if (tmp >= prime) tmp -= prime;
	if (tmp >= prime) tmp -= prime;

	int64_t res = multiplyReduce(a0, b0, prime) + double30Times(tmp, prime);
	if (res >= prime) res -= prime;

	return res;

#endif
}

int64_t inverse (int64_t a, int64_t prime) noexcept {
	return fastExp(a, prime - 2, prime);
}

int64_t fastExp (int64_t a, int64_t b, int64_t prime) noexcept {
	if (!b) return 1;
	int64_t tmp = fastExp(a, b >> 1, prime);
	tmp = multiply(tmp, tmp, prime);
	if (b & 1) tmp = multiply(tmp, a, prime);
	return tmp;
}

std::vector<std::vector<int64_t>> matrixMultiply (const std::vector<std::vector<int64_t>> &a,
		const std::vector<std::vector<int64_t>> &b, int64_t prime) {

	size_t n = a.size();
	std::vector<std::vector<int64_t>> res(n, std::vector<int64_t>(n));

	for (size_t i = 0; i < n; i++) {
		for (size_t k = 0; k < n; k++) {
			for (size_t j = 0; j < n; j++) {
				res[i][j] += multiply(a[i][k], b[k][j], prime);
				if (res[i][j] >= prime) res[i][j] -= prime;
			}
		}
	}

	return res;
}

std::vector<std::vector<int64_t>> matrixInverse (std::vector<std::vector<int64_t>> a, int64_t prime) {
	size_t n = a.size();
	std::vector<std::vector<int64_t>> res(n, std::vector<int64_t>(n));
	for (size_t i = 0; i < n; i++) res[i][i] = 1;

	auto swapRows = [&a, &res] (size_t i, size_t j) {
		std::swap(a[i], a[j]);
		std::swap(res[i], res[j]);
	};

	auto multiplyRow = [&a, &res, n, prime] (size_t i, int64_t c) {
		for (size_t j = 0; j < n; j++) {
			if (a[i][j]) a[i][j] = multiply(a[i][j], c, prime);
			if (res[i][j]) res[i][j] = multiply(res[i][j], c, prime);
		}
	};

	auto makeColumnUnit = [&a, &res, n, prime] (size_t i) { /* assumes a[i][i] = 1 */
		for (size_t k = 0; k < n; k++) { /* iterate over other rows */
			if (k == i || !a[k][i]) continue;
			for (size_t j = i + 1; j < n; j++) { /* iterate over columns in a */
				a[k][j] -= multiply(a[i][j], a[k][i], prime);
				if (a[k][j] < 0) a[k][j] += prime;
			}
			for (size_t j = 0; j < n; j++) { /* iterate over columns in res */
				res[k][j] -= multiply(res[i][j], a[k][i], prime);
				if (res[k][j] < 0) res[k][j] += prime;
			}
			a[k][i] = 0;
		}
	};

	for (size_t i = 0; i < n; i++) { /* make a 1 in row i, column i, assuming that the submatrix 0, …, i-1 is already unit */
		size_t row = n;
		for (size_t k = i; k < n; k++) {
			if (a[k][i]) {
				row = k;
				break;
			}
		}
		if (row == n) {
			throw std::runtime_error("matrixInverse: Matrix not invertible. Please report this error.");
		}
		swapRows(row, i);
		multiplyRow(i, inverse(a[i][i], prime));
		makeColumnUnit(i);
	}
	return res;
}

}
