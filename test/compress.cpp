#include <iostream>
#include <vector>

#include "template.hpp"
#include "compress.hpp"

signed main() {
    int n; std::cin >> n;
    std::vector<int> a(n);
    REP(i, n) std::cin >> a[i];
    Compress comp_a(ALL(a));
    debug(comp_a);
    REP(i, n) debug(comp_a(comp_a[i]));
    return 0;
}
