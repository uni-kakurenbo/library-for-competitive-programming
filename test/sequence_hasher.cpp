#include <bits/stdc++.h>
#include "template/standard.hpp"
#include "hash/sequence_hasher.hpp"

// https://atcoder.jp/contests/abc141/tasks/abc141_e

signed main() {
    int n; std::cin >> n;
    std::string s; std::cin >> s;

    lib::sequence_hasher hash(ALL(s));
    debug(lib::sequence_hasher<>::base);

    int ans = 0;
    REP(i, s.size()) REP(j, i, (int)s.size()) {
        while(i + ans < j and j + ans < n) {
            if(hash(i, i+ans+1) != hash(j, j+ans+1)) break;
            ++ans;
        }
    }

    print(ans);

    return 0;
}
