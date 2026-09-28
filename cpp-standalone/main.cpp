#include "../src/htc.h"

#include <cstring>
#include <iostream>
#include <stdexcept>

#if defined(HAVE_GMP)
#include <gmp.h>
#elif defined(HAVE_OPENSSL)
#include <openssl/bn.h>
#endif

namespace {

enum Verbosity {
	minimal,
	normal
};

#if defined(HAVE_GMP)

bool isProbablyPrimeGMP (int64_t candidate) {
	mpz_t z;
	mpz_init(z);
	mpz_import(z, 1, -1, sizeof(candidate), 0, 0, &candidate);
	int result = mpz_probab_prime_p(z, 40);
	mpz_clear(z);
	return result > 0;
}

#elif defined(HAVE_OPENSSL)

bool isProbablyPrimeOPENSSL (int64_t candidate) {
	unsigned char bytes[8];
	for (size_t i = 0; i < 8; i++) {
		bytes[i] = static_cast<unsigned char>((candidate >> (8 * i)) & 0xFF);
	}
	BIGNUM *bn = BN_lebin2bn(bytes, 8, nullptr);
	if (!bn) {
		throw std::runtime_error("OpenSSL: BN_lebin2bn failed");
	}
	BN_CTX *ctx = BN_CTX_new();
	if (!ctx) {
		BN_free(bn);
		throw std::runtime_error("OpenSSL: BN_CTX_new failed");
	}

#if OPENSSL_VERSION_NUMBER >= 0x30000000L
	int result = BN_check_prime(bn, ctx, nullptr);
#else
	int result = BN_is_prime_ex(bn, 40, ctx, nullptr);
#endif

	BN_CTX_free(ctx);
	BN_free(bn);
	if (result < 0) {
		throw std::runtime_error("OpenSSL prime check failed with error code < 0");
	}
	return result > 0;
}

#else

bool isProbablyPrimeError ([[maybe_unused]] int64_t candidate) {
	throw std::logic_error("Neither HAVE_GMP nor HAVE_OPENSSL defined but no prime given");
}

#endif

}

int main (int argc, char **argv) {
	std::optional<uint64_t> seed;
	std::optional<int64_t> prime;
	Verbosity verbosity = normal; 

	for (int i = 1; i < argc; i++) {
		if (std::strcmp(argv[i], "--seed") == 0) {
			if (! (i + 1 < argc)) {
				throw std::invalid_argument("When using --seed SEED, SEED cannot be omitted");
			}
			seed = std::stoull(argv[i + 1]);
			i++;
		} else if (std::strcmp(argv[i], "--prime") == 0) {
			if (! (i + 1 < argc)) {
				throw std::invalid_argument("When using --prime PRIME, PRIME cannot be omitted");
			}
			prime = std::stoll(argv[i + 1]);
			if (! ((1ll << 58) <= prime && prime < (1ll << 59))) {
				throw std::invalid_argument("When using --prime PRIME, PRIME must be between 2^58 and 2^59");
			}
			i++;
		} else if (std::strcmp(argv[i], "--minimal") == 0) {
			verbosity = minimal;
		} else if (std::strcmp(argv[i], "--help") == 0) {
			std::cout << "Usage: " << argv[0] << " [OPTION]..." << std::endl
					  << "Reads the input from stdin and writes to stdout" << std::endl << std::endl
					  << "Input format:" << std::endl
					  << "n m_1 m_2\nu_0 v_0\n…\nu_{m_1-1} v_{m_1-1}\ni_0 j_0\n…\ni_{m_2-1} j_{m_2-1}\n" << std::endl
					  << "where n is the number of nodes, m_1 is the number of bidirected edges," << std::endl
					  << "m_2 is the number of directed edges, the nodes are numbered 0, …, n-1," << std::endl
					  << "{u_0, v_0}, …, {u_{m_1-1},v_{m_1-1}} are the bidirected edges" << std::endl
					  << "and (i_0, j_0), …, (i_{m_2-1},j_{m_2-1}) are the directed edges" << std::endl << std::endl
					  << "Options:" << std::endl
					  << "   --seed SEED        Seed the random with SEED. Generates a random seed by default." << std::endl
					  << "   --prime PRIME      Use PRIME as the prime. PRIME should be between 2^58 and 2^59." << std::endl
					  << "                      Selects a random prime in this range by default." << std::endl
					  << "   --minimal          Only output the identified directed edges (i, j)" << std::endl
					  << "                      in the format i j, one edge per line" << std::endl
					  << "   --help             Display this help and exit" << std::endl;
			return 0;
		} else {
			throw std::invalid_argument(std::string("Unknown command line option: ") + argv[i]);
		}
	}

	size_t n, m1, m2;
	if (! ((std::cin >> n) && n > 0) || ! (std::cin >> m1) || ! (std::cin >> m2)) {
		throw std::invalid_argument("The first line of the input should contain three integers n m_1 m_2, the number of nodes, the number of bidirected edges and the number of directed edges");
	}
	std::vector<std::pair<size_t, size_t>> bidirected(m1), directed(m2);
	for (size_t i = 0; i < m1; i++) {
		size_t u, v;
		if (! ((std::cin >> u >> v) && u != v && u < n && v < n)) {
			throw std::invalid_argument("Lines 2 to 2+m_1-1 of the input should each contain two integers u and v with 0<=u,v<n and u!=v, indicating that there is a bidirected edge between u and v");
		}
		if (u > v) std::swap(u, v);
		bidirected[i] = {u, v};
	}
	for (size_t i = 0; i < m2; i++) {
		size_t u, v;
		if (! ((std::cin >> u >> v) && u != v && u < n && v < n)) {
			throw std::invalid_argument("Lines 2+m_1 to 2+m_1+m_2-1 of the input should each contain two integers i and j with 0<=i,j<n and i!=j, indicating that there is a directed edge from i to j");
		}
		directed[i] = {u, v};
	}

#if defined(HAVE_OPENSSL)
	auto identificationResult = fasterhtc::htc(n, bidirected, directed, seed, prime, isProbablyPrimeOPENSSL);
#elif defined(HAVE_GMP)
	auto identificationResult = fasterhtc::htc(n, bidirected, directed, seed, prime, isProbablyPrimeGMP);
#else
	auto identificationResult = fasterhtc::htc(n, bidirected, directed, seed, prime, isProbablyPrimeError);
#endif

	if (verbosity != minimal) {
		std::cout << "Identification completed with seed = " << identificationResult.seed << " and prime = " << identificationResult.prime << std::endl << "Identified the following " << identificationResult.identification.size() << " edges:" << std::endl;
	}
	for (auto [i, j] : identificationResult.identification) {
		if (verbosity == minimal) {
			std::cout << i << " " << j << std::endl;
		} else {
			std::cout << i << " -> " << j << std::endl;
		}
	}
}
