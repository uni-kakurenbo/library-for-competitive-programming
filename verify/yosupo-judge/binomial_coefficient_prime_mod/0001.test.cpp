/*
 * @uni_kakurenbo
 * https://github.com/uni-kakurenbo/competitive-programming-workspace
 *
 * CC0 1.0  http://creativecommons.org/publicdomain/zero/1.0/deed.ja
 */
/* #language C++ 20 GCC */

#define PROBLEM "https://judge.yosupo.jp/problem/binomial_coefficient_prime_mod"

#include "sneaky/enforce_int128_enable.hpp"

#include <iostream>

#include "snippet/aliases.hpp"
#include "snippet/fast_io.hpp"
#include "snippet/iterations.hpp"
#include "adaptor/io.hpp"
#include "snippet/aliases.hpp"
#include "numeric/binomial_coefficient.hpp"

signed main() {
    int t; std::cin >> t;
    lib::i32 m; std::cin >> m;

    lib::binomial_coefficient_prime_power_mod<lib::i32> binom(m, 1, 10'000'000);

    REP(t) {
        lib::i32 n, k; std::cin >> n >> k;
        print(binom.comb(n, k));
    }
}
