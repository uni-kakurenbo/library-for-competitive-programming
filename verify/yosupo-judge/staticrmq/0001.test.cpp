/*
 * @uni_kakurenbo
 * https://github.com/uni-kakurenbo/competitive-programming-workspace
 *
 * CC0 1.0  http://creativecommons.org/publicdomain/zero/1.0/deed.ja
 */
/* #language C++ GCC */

#define PROBLEM "https://judge.yosupo.jp/problem/staticrmq"

#include <iostream>
#include "snippet/aliases.hpp"
#include "snippet/fast_io.hpp"
#include "snippet/iterations.hpp"
#include "adaptor/io.hpp"
#include "adaptor/valarray.hpp"
#include "data_structure/disjoint_sparse_table.hpp"
#include "action/range_min.hpp"


signed main() {
    int n, q; std::cin >> n >> q;
    lib::valarray<int> a(n); input >> a;

    lib::disjoint_sparse_table<lib::actions::range_min<int>> min(a);

    REP(q) {
        int l, r; std::cin >> l >> r;
        print(min(l, r).fold());
    }
}
