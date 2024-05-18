#pragma once


#include <vector>
#include <ranges>
#include <stack>
#include <bit>


#include "snippet/aliases.hpp"

#include "internal/dev_env.hpp"

#include "data_structure/disjoint_sparse_table.hpp"
#include "algebraic/helper.hpp"
#include "numeric/limits.hpp"


namespace uni {


template<class Graph>
struct compressed_tree {
    using size_type = Graph::size_type;

  private:
    std::vector<size_type> _in, _depth, _h;
    std::vector<std::vector<size_type>> _table;


    void _dfs(const Graph& tree, const size_type v, const size_type p = -1, const size_type depth = 0) noexcept(NO_EXCEPT) {
        this->_depth[v] = depth;

        this->_in[v] = this->_table[0].size();
        this->_table[0].push_back(v);

        ITR(e, tree[v]) {
            if(e == p) continue;
            this->_dfs(tree, e, v, depth + 1);
            this->_table[0].push_back(v);
        }

        this->_table[0].push_back(v);
    }

  public:
    compressed_tree(const Graph& tree, const size_type root = 0) noexcept(NO_EXCEPT)
      : _in(tree.size()), _depth(tree.size()), _table(1)
    {
        this->_dfs(tree, root);

        const auto size = this->_table[0].size();

        this->_h.assign(size + 1, 0);
        FOR(i, 2, size) this->_h[i] = this->_h[i >> 1] + 1;
        debug(this->_h);

        this->_table.resize(this->_h[size] + 1);
        FOR(i, 1, this->_h[size]) {
            this->_table[i].assign(size - (1 << i) + 1, size);
        }

        size_type d = 1;
        REP(i, this->_h[size]) {
            FOR(j, size - (d << 1)) {
                this->_table[i + 1][j] =
                    std::min(
                        this->_table[i][j],
                        this->_table[i][j + d],
                        [&](const auto a, const auto b) { return this->_depth[a] < this->_depth[b]; }
                    );
            }
            d <<= 1;
        }

        debug(this->_table);
    }


    template<std::ranges::input_range R, class Res>
        requires std::convertible_to<std::ranges::range_value_t<R>, size_type>
    std::optional<size_type> build(R&& ranges, Res *const res) const noexcept(NO_EXCEPT) {
        if(std::ranges::size(ranges) == 0) return {};

        vector<size_type> vs(ALL(ranges));
        std::ranges::sort(vs, [&](const auto i, const auto j) noexcept(NO_EXCEPT) { return this->_in[i] < this->_in[j]; });

        std::vector<size_type> stack;
        stack.push_back(vs[0]);

        res->clear();

        REP(i, std::ranges::ssize(vs) - 1) {
            const auto w = this->lca(vs[i], vs[i + 1]);

            if(w != vs[i]) {
                auto last = stack.back(); stack.pop_back();

                while(!stack.empty() && this->_depth[w] <= this->_depth[stack.back()]) {
                    res->operator[](last) = stack.back();
                    last = stack.back(); stack.pop_back();
                }

                if(stack.empty() || stack.back() != w) {
                    stack.push_back(w);
                    vs.push_back(w);
                }
                res->operator[](last) = w;
            }

            stack.push_back(vs[i + 1]);
        }

        REP(i, std::ranges::ssize(stack) - 1) {
            res->operator[](stack[i + 1]) = stack[i];
        }

        return stack.front();
    }

    auto lca(const size_type v0, const size_type v1) const noexcept(NO_EXCEPT) {
        assert(0 <= v0 && v0 < std::ranges::ssize(this->_in));
        assert(0 <= v1 && v1 < std::ranges::ssize(this->_in));

        auto p0 = this->_in[v0], p1 = this->_in[v1];
        if(p0 > p1) std::swap(p0, p1);

        const auto d = this->_h[p1 - p0 + 1];
        debug(d, p0, p1);

        return
            std::min(
                this->_table[d][p0],
                this->_table[d][p1 - (1 << d) + 1],
                [&](const auto a, const auto b) { return this->_depth[a] < this->_depth[b]; }
            );
    }

    auto depth(const size_type v) const noexcept(NO_EXCEPT) {
        assert(0 <= v && v <= std::ranges::ssize(this->_depth));
        return this->_depth[v];
    }
};


} // namespace uni
