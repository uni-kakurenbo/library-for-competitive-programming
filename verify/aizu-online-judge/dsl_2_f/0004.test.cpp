/*
 * @uni_kakurenbo
 * https://github.com/uni-kakurenbo/competitive-programming-workspace
 *
 * CC0 1.0  http://creativecommons.org/publicdomain/zero/1.0/deed.ja
 */
/* #language C++ 20 GCC */

#define PROBLEM "https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=DSL_2_F"

#include <iostream>
#include "snippet/aliases.hpp"
#include "snippet/fast_io.hpp"
#include "snippet/iterations.hpp"
#include "adaptor/io.hpp"
#include "data_structure/dynamic_sequence.hpp"
#include "action/range_min.hpp"
#include "action/helpers.hpp"

signed main() {
    int n, q; input >> n >> q;

    lib::dynamic_sequence<lib::actions::make_full_t<lib::actions::range_min<int>>> data(n);

    REP(q) {
        int t, l, r; input >> t >> l >> r; ++r;
        if(t == 0) {
            int x; input >> x;
            data(l, r) = x;
        }
        if(t == 1) {
            print(data(l, r).fold());
        }
    }
}
