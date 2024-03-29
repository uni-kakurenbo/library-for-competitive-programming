/*
 * @uni_kakurenbo
 * https://github.com/uni-kakurenbo/competitive-programming-workspace
 *
 * CC0 1.0  http://creativecommons.org/publicdomain/zero/1.0/deed.ja
 */
/* #language C++ 20 GCC */

#define PROBLEM "https://judge.yosupo.jp/problem/longest_increasing_subsequence"

#include <iostream>
#include "snippet/fast_io.hpp"
#include "adaptor/io.hpp"
#include "adaptor/vector.hpp"
#include "iterable/longest_increasing_subsequence.hpp"

signed main() {
    int n; std::cin >> n;
    lib::vector<int> a(n); input >> a;
    lib::lis<true,int> lis(ALL(a));

    print(lis.size(), lis.indices);
}
