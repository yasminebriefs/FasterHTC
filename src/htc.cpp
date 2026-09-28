#include "htc.h"

#include <chrono>
#include <deque>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>

namespace fasterhtc {

std::vector<std::vector<int64_t>> computeSigmaInTermsOfLambdaAndOmega (std::vector<std::vector<int64_t>> &lambda,
		const std::vector<std::vector<int64_t>> &omega, int64_t prime) {

	/* Works in O(n^3) */

	size_t n = lambda.size();

	/* Temporarily compute I - Λ */
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (lambda[i][j]) lambda[i][j] = prime - lambda[i][j]; /* In range [1, prime-1], so no modulo needed */
		}
	}
	for (size_t i = 0; i < n; i++) {
		lambda[i][i] = 1;
	}

	auto lambdaInverse = matrixInverse(lambda, prime); /* Compute (I - Λ)^(-1) */

	/* Restore Λ */
	for (size_t i = 0; i < n; i++) {
		lambda[i][i] = 0;
	}
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (lambda[i][j]) lambda[i][j] = prime - lambda[i][j];
		}
	}

	auto tmp = matrixMultiply(omega, lambdaInverse, prime);

	/* Transpose lambdaInverse */
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i+1; j < n; j++) {
			std::swap(lambdaInverse[i][j], lambdaInverse[j][i]);
		}
	}

	return matrixMultiply(lambdaInverse, tmp, prime);
}

void computeHalfTreks (const std::vector<std::pair<size_t, size_t>> &bidirected,
		const std::vector<std::vector<size_t>> &parents,
		std::vector<std::vector<bool>> &halfTrek) {

	/* For each node, compute all nodes reachable by half-treks in O(n^2 + m2 * n)) ⊂ O(n^3) */

	size_t n = halfTrek.size();

	std::vector<std::vector<size_t>> adj(n);
	for (size_t j = 0; j < n; j++) {
		for (size_t i : parents[j]) {
			if (i < n) {
				adj[i].push_back(j);
			}
		}
	}

	std::function<void (size_t, std::vector<bool>&)> computeHalfTreksDFS = [&adj, &computeHalfTreksDFS]
		(size_t u, std::vector<bool> &halfTrekV) {

		if (halfTrekV[u]) return;
		halfTrekV[u] = true;
		for (size_t v : adj[u]) computeHalfTreksDFS(v, halfTrekV);
	};

	for (size_t u = 0; u < n; u++) {
		computeHalfTreksDFS(u, halfTrek[u]);
	}
	for (auto [u, v] : bidirected) {
		computeHalfTreksDFS(v, halfTrek[u]);
		computeHalfTreksDFS(u, halfTrek[v]);
	}
}

bool addRow (int64_t prime, std::vector<int64_t> &&row,
		std::vector<std::optional<std::vector<int64_t>>> &rows) {

	for (size_t i = 0; i < row.size(); i++) {
		/* Invariant: row[0 ...  i-1] = 0 */

		if (row[i] == 0) continue;
		if (rows[i]) {
			/* Make the entry in this column zero */
			for (size_t j = row.size() - 1; j > i; j--) {
				row[j] += prime - multiply(row[i], (*rows[i])[j], prime);
				if (row[j] >= prime) row[j] -= prime;
			}
			row[i] = 0;
			continue;
		}
		/* Multiply the row to set the leading entry to 1 */
		int64_t inv = inverse(row[i], prime);
		row[i] = 1;
		for (size_t j = i + 1; j < row.size(); j++) {
			row[j] = multiply(row[j], inv, prime);
		}
		rows[i].emplace(std::move(row));
		return true;
	}
	return false;
}

std::vector<int64_t> computeRow (size_t y, size_t v, size_t numBasic, int64_t prime,
		const std::vector<std::vector<int64_t>> &lambda,
		const std::vector<std::vector<int64_t>> &sigma,
		const std::vector<std::vector<size_t>> &parents) {

	/* This function is only called for basic nodes v */
	std::vector<int64_t> row(parents[v].size());
	for (size_t i = 0; i < parents[v].size(); i++) {
		size_t w = parents[v][i];
		row[i] = sigma[y][w];
		if (y < numBasic) { /* Non-basic nodes don't have parents */
			for (size_t h : parents[y]) {
				row[i] += prime - multiply(sigma[h][w], lambda[h][y], prime);
				if (row[i] >= prime) row[i] -= prime;
			}
		}
	}

	return row;
}

void identifyComponent (size_t numBasic, const std::vector<size_t> &translation, const std::shared_ptr<RandomTool> &randomTool,
		const std::vector<std::pair<size_t, size_t>> &bidirected,
		const std::vector<std::vector<size_t>> &parents,
		IdentificationResult &identificationResult) {

	size_t totalNum = translation.size();
	std::vector<std::vector<int64_t>> lambda(totalNum, std::vector<int64_t>(totalNum)),
										omega(totalNum, std::vector<int64_t>(totalNum));

	for (size_t j = 0; j < numBasic; j++) {
		for (size_t i : parents[j]) {
			lambda[i][j] = randomTool->getRandomSmallerPrime();
		}
	}
	for (auto [u, v] : bidirected) {
		omega[u][v] = randomTool->getRandomSmallerPrime();
		omega[v][u] = omega[u][v];
	}
	for (size_t i = 0; i < totalNum; i++) {
		omega[i][i] = randomTool->getRandomSmallerPrime();
	}

	auto sigma = computeSigmaInTermsOfLambdaAndOmega(lambda, omega, randomTool->getPrime());

	std::vector<std::vector<bool>> halfTrek(numBasic, std::vector<bool>(numBasic));
	computeHalfTreks(bidirected, parents, halfTrek); /* The extra nodes are not reachable by half treks */

	/* Now perform the identification */
	std::vector<size_t> numColumns(numBasic);
	std::vector<std::vector<std::optional<std::vector<int64_t>>>> rows(numBasic);
	for (size_t i = 0; i < numBasic; i++) {
		rows[i].assign(parents[i].size(), std::nullopt);
	}

	for (size_t y = 0; y < totalNum; y++) {
		for (size_t v = 0; v < numBasic; v++) {
			if (numColumns[v] == parents[v].size()) continue; /* v is already identified */
			if (y >= numBasic || !halfTrek[v][y]) {
				/* y ∉ htr(v) implies that y is not v nor a sibling of v */
				if (addRow(identificationResult.prime, computeRow(y, v, numBasic, identificationResult.prime, lambda, sigma, parents), rows[v])) {
					numColumns[v]++;
				}
			}
		}
	}

	std::deque<size_t> queue;
	for (size_t i = 0; i < numBasic; i++) {
		if (numColumns[i] == parents[i].size()) {
			queue.push_back(i);
		}
	}


	while (!queue.empty()) {
		size_t y = queue.front();
		queue.pop_front();

		for (size_t p : parents[y]){
			identificationResult.addIdentifiedEdgeTranslated(p, y, translation);
		}

		if (y >= numBasic) continue; /* y's row was already added for all v */
		for (size_t v = 0; v < numBasic; v++) {
			if (numColumns[v] == parents[v].size()) continue; /* v is already identified */
			if (omega[v][y]) continue; /* y can't be v or a sibling of v */
			if (halfTrek[v][y]) { /* The row was not added before */
				if (addRow(identificationResult.prime, computeRow(y, v, numBasic, identificationResult.prime, lambda, sigma, parents), rows[v])) {
					numColumns[v]++;
					if (numColumns[v] == parents[v].size()) {
						queue.push_back(v);
					}
				}
			}
		}
	}
}

size_t unionFindRoot (std::vector<int> &unionFind, size_t u) {
	if (unionFind[u] < 0) return u;
	return unionFind[u] = unionFindRoot(unionFind, unionFind[u]);
}

void unionFindMerge (std::vector<int> &unionFind, size_t u, size_t v) {
	u = unionFindRoot(unionFind, u);
	v = unionFindRoot(unionFind, v);
	if (u == v) return;
	if (unionFind[u] > unionFind[v]) std::swap(u, v);
	unionFind[u] += unionFind[v];
	unionFind[v] = u;
}

void mergeStronglyConnectedComponents (size_t n, std::vector<int> &unionFind,
		const std::vector<std::pair<size_t, size_t>> &directed) {

	std::vector<std::vector<size_t>> graph(n);

	for (auto [u, v] : directed) {
		graph[u].push_back(v);
	}

	std::vector<int> low(n), num(n);
	std::vector<size_t> stack;
	std::vector<bool> active(n);
	size_t curTime = 0;

	std::function<void (size_t)> dfsSCC = [&low, &num, &stack, &active, &curTime, &graph, &unionFind, &dfsSCC] (size_t u) {
		low[u] = num[u] = ++curTime;
		stack.push_back(u);
		active[u] = true;
		for (size_t v : graph[u]) {
			if (!num[v]) dfsSCC(v);
			if (active[v]) low[u] = std::min(low[u], low[v]);
		}
		if (low[u] == num[u]) {
			size_t v;
			do {
				v = stack.back();
				stack.pop_back();
				active[v] = false;
				unionFindMerge(unionFind, u, v);
			} while (u != v);
		}
	};

	for (size_t u = 0; u < n; u++) {
		if (!num[u]) dfsSCC(u);
	}
}


void solveByTianDecomposition (size_t n, IdentificationResult &identificationResult, const std::shared_ptr<RandomTool> &randomTool,
		const std::vector<std::pair<size_t, size_t>> &bidirected,
		const std::vector<std::pair<size_t, size_t>> &directed) {

	std::vector<int> unionFind(n, -1); /* -x means points to itself, x nodes in subtree, +y means points to y */

	mergeStronglyConnectedComponents(n, unionFind, directed);

	for (auto [u, v] : bidirected) {
		unionFindMerge(unionFind, u, v);
	}

	std::vector<std::vector<size_t>> parents(n);

	for (auto [u, v] : directed) {
		parents[v].push_back(u);
	}

	std::vector<size_t> indexInComponent(n, n); /* n means not in the component */

	for (size_t i = 0; i < n; i++) { /* Iterate over the components */
		if (unionFindRoot(unionFind, i) != i) continue;
		std::vector<size_t> component;
		for (size_t j = 0; j < n; j++) {
			if (unionFindRoot(unionFind, j) == i) {
				component.push_back(j);
				indexInComponent[j] = 0;
			}
		}
		size_t numBasic = component.size();
		bool hasEdges = false;
		for (size_t j = 0; j < numBasic; j++) {
			size_t v = component[j];
			for (size_t w : parents[v]) {
				hasEdges = true;
				if (indexInComponent[w] == n) {
					component.push_back(w);
					indexInComponent[w] = 0;
				}
			}
		}
		if (!hasEdges){
			for (size_t j : component) {
				indexInComponent[j] = n;
			}
			continue; /* No directed edges to identify in this component */
		}
		for (size_t j = 0; j < component.size(); j++) {
			indexInComponent[component[j]] = j;
		}
		std::vector<std::pair<size_t, size_t>> tmpBidirected;
		std::vector<std::vector<size_t>> tmpParents(numBasic);
		for (auto [u, v] : bidirected) {
			if (unionFindRoot(unionFind, u) == i && unionFindRoot(unionFind, v) == i) {
				tmpBidirected.emplace_back(indexInComponent[u], indexInComponent[v]);
			}
		}
		for (size_t j = 0; j < numBasic; j++) {
			size_t v = component[j];
			for (size_t w : parents[v]) {
				tmpParents[j].push_back(indexInComponent[w]);
			}
		}
		identifyComponent(numBasic, component, randomTool, tmpBidirected, tmpParents, identificationResult);

		for (size_t j : component) {
			indexInComponent[j] = n;
		}
	}
}

uint64_t generateSeed () {
	std::random_device rd;
	uint64_t rd_seed = (static_cast<uint64_t>(rd()) << 32) | rd(); /* May be deterministic on some platforms */
	uint64_t time_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	return rd_seed ^ time_seed; /* So xor with the current time to get something "non-deterministic" on all platforms */
}

IdentificationResult htc (size_t n, const std::vector<std::pair<size_t, size_t>> &bidirected,
		const std::vector<std::pair<size_t, size_t>> &directed, std::optional<uint64_t> seed_opt,
		std::optional<int64_t> prime_opt, bool (*isProbablyPrime)(int64_t)) {

	uint64_t seed = seed_opt.has_value() ? seed_opt.value() : generateSeed();
	auto randomTool = std::make_shared<RandomTool>(seed, prime_opt, isProbablyPrime);
	int64_t prime = randomTool->getPrime();

	/* Validate input */
	if (!((1ll << 58) <= prime && prime < (1ll << 59))) {
		throw std::invalid_argument("The prime should be between 2^58 and 2^59");
	}

	for (auto [i, j] : directed) {
		if (!(i < n && j < n && i != j)) {
			throw std::invalid_argument("Invalid input: for directed edges i->j, it must hold that 0 <= i, j < n and i != j");
		}
	}

	for (auto [u, v] : bidirected) {
		if (!(u < n && v < n && u != v)) {
			throw std::invalid_argument("Invalid input: for bidirected edges u<->v, it must hold that 0 <= u, v < n and u != v");
		}
	}

	IdentificationResult identificationResult(seed, prime);

	solveByTianDecomposition(n, identificationResult, randomTool, bidirected, directed);

	return identificationResult;
}

}
