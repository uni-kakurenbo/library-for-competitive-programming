#include <iostream>

#include "snippet/iterations.hpp"

#include "data_structure/range_operation/range_affine_range_minmax.hpp"

#include "template.hpp"

signed main() {
    int n, q; std::cin >> n >> q;

    Lib::RangeAffineRangeMinMax<ll> data;
    data.resize(n);
    debug(data.size());

    REP(q) {
        int t; std::cin >> t;
        if(t == 0) {
            int l, r, x; std::cin >> l >> r >> x; ++r;
            data.add(l, r, x);
        }
        if(t == 1) {
            int l, r; std::cin >> l >> r; ++r;
            std::cout << data.min(l, r) << "\n";
        }
        // debug(data);
    }
}
