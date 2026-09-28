FILES := main.cpp ../src/htc.cpp ../src/algebra.cpp ../src/random.cpp

OPENSSL := -lcrypto -DHAVE_OPENSSL
GMP := -lgmp -DHAVE_GMP

default: with-openssl

with-openssl:
	g++ -std=c++17 -O2 $(FILES) $(OPENSSL) -o FasterHTC

with-gmp:
	g++ -std=c++17 -O2 $(FILES) $(GMP) -o FasterHTC

no-prime-check:
	g++ -std=c++17 -O2 $(FILES) -o FasterHTC

debug:
	g++ -std=c++17 -g -Wall -Wextra -Wshadow -fsanitize=undefined,address -D_GLIBCXX_DEBUG $(FILES) $(OPENSSL) -o FasterHTC

clean:
	rm FasterHTC
