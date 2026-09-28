# FasterHTC
This is an implementation of the algorithm from the paper 'A Faster Algorithm for the Half-Trek Criterion in Structural Causal Models' (2026) by Briefs and Bläser.
Given a structural causal model (SCM), it determines which parameters are identifiable by the half-trek criterion.
The algorithm is randomized and its running time is $O(n^3+n^2m_2)$ where $n$ is the number of nodes and $m_2$ is the number of directed edges.
This simplifies to $O(n^4)$ for sparse and $O(n^5)$ for dense graphs.

## The C++ program
The C++ program can be compiled by running `make` in cpp-standalone/.
This requires OpenSSL to be installed (install by `apt install libssl-dev` in Ubuntu for example).
Alternatively, using `make with-gmp`, the program can be compiled with GMP instead of OpenSSL.
By `make no-prime-check`, the program can be compiled without OpenSSL and GMP, but then, a 59-bit prime has to be given manually every call.
Otherwise, the program can generate a random prime using GMP or OpenSSL.

The executable is called `FasterHTC`.
It reads the input from stdin and writes to stdout.

**Input format**

$n$ $m_1$ $m_2$  
$u_0$ $v_0$  
$\dots$  
$u_{m_1-1}$ $v_{m_1-1}$  
$i_0$ $j_0$  
$\dots$  
$i_{m_2-1}$ $v_{m_2-1}$  

where $n$ is the number of nodes, $m_1$/$`m_2`$ is the number of bidirected/directed edges,
the nodes are numbered $0,\dots, n-1$,
$`\{u_0, v_0\}, \dots, \{u_{m_1-1},v_{m_1-1}\}`$ are the bidirected edges, and
$(i_0, j_0), \dots, (i_{m_2-1},j_{m_2-1})$ are the directed edges

**Options**
```
   --seed SEED        Seed the random with SEED. Generates a random seed by default.
   --prime PRIME      Use PRIME as the prime. PRIME should be between 2^58 and 2^59.
                      Selects a random prime in this range by default.
   --minimal          Only output the identified directed edges (i, j)
                      in the format i j, one edge per line
   --help             Display this help and exit
```

**Example**
```
> ./FasterHTC < ../test-files/handmade/simple_4_nodes_in.txt
Identification completed with seed = 15788885456488170130 and prime = 477851502274213721
Identified the following 3 edges:
0 -> 2
1 -> 3
2 -> 3
```

## General notes
The algorithm is randomized.
The error probability is below $3.52e-7$ for $n=1000$ and below $3.69e-4$ for $n=10\,000$.
These estimates are very conservative worst-case estimates, and the algorithm can be run multiple times with different seeds (and primes) to decrease the error probability exponentially.

If no seed is given, a seed is generated randomly.
If the prime is not given, using the same seed (on the same system) twice will generate the same prime twice.
Whether or not a prime is given doesn't influence the other random computations.

The program is based on [FastTreeID](https://github.com/yasminebriefs/FastTreeID), an implementation of a complete algorithm for generic identification in tree-shaped SCMs.

## References
Yasmine Briefs and Markus Bläser. A faster algorithm for the half-trek criterion in structural causal models. Advances in Neural Information Processing Systems 30: Annual Conference on Neural Information Processing Systems 2026, NeurIPS 2026, December 6-12, 2026, Sydney, Australia (To appear)

Yasmine Briefs and Markus Bläser. Faster generic identification in tree-shaped structural causal models. Advances in Neural Information Processing Systems 29: Annual Conference on Neural Information Processing Systems 2025, NeurIPS 2025, December 2-7, 2025, San Diego, USA [https://openreview.net/forum?id=8PHOPPH35D](https://openreview.net/forum?id=8PHOPPH35D)

Rina Foygel, Jan Draisma, and Mathias Drton. Half-trek criterion for generic identifiability of linear structural equation models. The Annals of Statistics, pages 1682--1713, 2012 [https://doi.org/10.1214/12-AOS1012](https://doi.org/10.1214/12-AOS1012)
