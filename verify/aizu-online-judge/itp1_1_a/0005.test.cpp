/*
 * @uni_kakurenbo
 * https://github.com/uni-kakurenbo/competitive-programming-workspace
 *
 * CC0 1.0  http://creativecommons.org/publicdomain/zero/1.0/deed.ja
 */
/* #language C++ GCC */

#define PROBLEM "https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=ITP1_1_A"

#include "sneaky/enforce_int128_enable.hpp"

#include <iostream>
#include "snippet/aliases.hpp"
#include "snippet/fast_io.hpp"
#include "snippet/iterations.hpp"
#include "adaptor/io.hpp"
#include "adaptor/set.hpp"
#include "data_structure/dynamic_set.hpp"
#include "data_structure/treap.hpp"
#include "action/null.hpp"
#include "random/engine.hpp"
#include "utility/timer.hpp"

signed main() {
    print("Hello World");

    lib::dynamic_set<lib::i64> data;
    lib::multiset<lib::i64> corr;

    debug(data);

    lib::timer timer(10000);

    while(not timer.expired()) {
        int t = lib::randi64() % 2;
        lib::i64 v = (lib::randi64() % 2) || data.empty() ? lib::randi64() : data[lib::randi64() % data.size()].val();
        debug(t, v);

        if(t == 0) {
            corr.insert(v);
            data.insert(v);
        }

        if(t == 1) {
            corr.remove(v);
            data.erase(v);
        }

        debug(corr, data);

        bool ok = true;
        {
            int i = 0;
            ITR(v, corr) {
                ok &= data[i++].val() == v;
            }
        }
        if(!ok) {
            assert(false);
        }
    }
}
