/*
 * @uni_kakurenbo
 * https://github.com/uni-kakurenbo/competitive-programming-workspace
 *
 * CC0 1.0  http://creativecommons.org/publicdomain/zero/1.0/deed.ja
 */
/* #language C++ GCC */
/* #region template */
#include "template.hpp"
/* #endregion */

#include "data_structure/fenwick_tree.hpp"
#include "data_structure/range_action/range_sum.hpp"


signed main() {
    int n, q; cin >> n >> q;
    vector<ll> a(n); input >> a;

    lib::fenwick_tree<lib::actions::range_sum<ll>> data(ALL(a));
    debug(data);

    REP(q) {
        int t; cin >> t;
        if(t == 0) {
            int p, x; cin >> p >> x;
            data[p] <<= x;
        }
        if(t == 1) {
            int l, r; cin >> l >> r;
            print(data(l, r).fold());
        }
        debug(data);
    }

    return 0;
}
