#pragma once


#include <cassert>
#include <algorithm>
#include <vector>
#include <type_traits>
#include <concepts>
#include <ranges>
#include <bit>


#include "snippet/iterations.hpp"

#include "internal/dev_env.hpp"
#include "internal/types.hpp"
#include "internal/iterator.hpp"
#include "internal/point_reference.hpp"
#include "internal/range_reference.hpp"
#include "internal/unconstructible.hpp"

#include "numeric/bit.hpp"
#include "algebraic/internal/concepts.hpp"
#include "action/base.hpp"


namespace lib {

namespace internal {

namespace segment_tree_impl {


// Thanks to: atcoder::segtree
template<algebraic::internal::monoid Monoid>
struct core {
    using size_type = internal::size_t;

    using operand = Monoid;

  protected:
    size_type _n = 0, _size = 0, _depth = 0;
    std::valarray<operand> _data;

    inline void pull(const size_type k) noexcept(NO_EXCEPT) {
        this->_data[k] = this->_data[k << 1] + this->_data[k << 1 | 1];
    }

  public:
    core() noexcept = default;

    explicit core(const size_type n) noexcept(NO_EXCEPT)
      : _n(n), _size(std::bit_ceil(lib::to_unsigned(n))), _depth(std::countr_zero(lib::to_unsigned(this->_size))),
        _data(this->_size << 1)
    {}


    inline void add(size_type p, const operand& x) noexcept(NO_EXCEPT) {
        this->set(p, this->_data[p + this->_size] + x);
    }

    inline void set(size_type p, const operand& x) noexcept(NO_EXCEPT) {
        p += this->_size;
        this->_data[p] = x;
        FOR(i, 1, this->_depth) this->pull(p >> i);
    }

    inline operand get(size_type p) const noexcept(NO_EXCEPT) {
        return this->_data[p + this->_size];
    }

    inline operand fold(size_type l, size_type r) const noexcept(NO_EXCEPT) {
        operand sml, smr;
        l += this->_size;
        r += this->_size;

        while(l < r) {
            if(l & 1) sml = sml + this->_data[l++];
            if(r & 1) smr = this->_data[--r] + smr;
            l >>= 1;
            r >>= 1;
        }
        return sml + smr;
    }

    inline operand fold_all() const noexcept(NO_EXCEPT) { return this->_data[1]; }
};



} // namespace segment_tree_impl

} // namespace internal


template<class T>
struct segment_tree : internal::unconstructible {};


template<algebraic::internal::monoid Monoid>
struct segment_tree<Monoid> : private internal::segment_tree_impl::core<Monoid> {
  private:
    using core = typename internal::segment_tree_impl::core<Monoid>;

  public:
    using value_type = Monoid;
    using size_type = typename core::size_type;

  private:
    inline size_type _positivize_index(const size_type p) const noexcept(NO_EXCEPT) {
        return p < 0 ? this->_n + p : p;
    }

  public:
    segment_tree() noexcept(NO_EXCEPT) : core() {};
    explicit segment_tree(const size_type n, const value_type& v = {}) noexcept(NO_EXCEPT) : core(n) { this->fill(v); }

    template<std::convertible_to<value_type> T>
    segment_tree(const std::initializer_list<T>& init_list) noexcept(NO_EXCEPT) : segment_tree(ALL(init_list)) {}

    template<std::input_iterator I, std::sized_sentinel_for<I> S>
    segment_tree(I first, S last) noexcept(NO_EXCEPT)
      : segment_tree(static_cast<size_type>(std::ranges::distance(first, last)))
    { this->assign(first, last); }

    template<std::ranges::input_range R>
    explicit segment_tree(R&& range) noexcept(NO_EXCEPT) : segment_tree(ALL(range)) {}


    inline size_type size() const noexcept(NO_EXCEPT) { return this->_n; }
    inline size_type allocated() const noexcept(NO_EXCEPT) { return this->_data.size(); }
    inline size_type depth() const noexcept(NO_EXCEPT) { return this->_depth; }


    template<std::convertible_to<value_type> T>
    inline auto& assign(const std::initializer_list<T>& init_list) noexcept(NO_EXCEPT) { return this->assign(ALL(init_list)); }

    template<std::input_iterator I, std::sentinel_for<I> S>
    inline auto& assign(I first, S last) noexcept(NO_EXCEPT) {
        if constexpr(std::sized_sentinel_for<S, I>) {
            assert(std::ranges::distance(first, last) == this->_n);
        }
        {
            size_type p = 0;
            for(auto itr=first; itr!=last; ++itr, ++p) this->_data[this->_size + p] = static_cast<value_type>(*itr);
        }
        REPD(p, 1, this->_size) this->pull(p);
        return *this;
    }

    template<std::ranges::input_range R>
    inline auto& assign(R&& range) noexcept(NO_EXCEPT) { return this->assign(ALL(range)); }

    inline auto& fill(const value_type& v = {}) noexcept(NO_EXCEPT) {
        REP(p, this->_n) this->_data[this->_size + p] = v;
        REPD(p, 1, this->_size) this->pull(p);
        return *this;
    }

    bool empty() const noexcept(NO_EXCEPT) { return this->_n == 0; }

    struct point_reference : internal::point_reference<segment_tree> {
        point_reference(segment_tree *const super, const size_type p) noexcept(NO_EXCEPT)
          : internal::point_reference<segment_tree>(super, super->_positivize_index(p))
        {
            assert(0 <= this->_pos && this->_pos < this->_super->_n);
        }

        operator value_type() const noexcept(NO_EXCEPT) { return this->_super->get(this->_pos); }
        value_type val() const noexcept(NO_EXCEPT) { return this->_super->get(this->_pos); }

        inline point_reference& set(const value_type& v) noexcept(NO_EXCEPT) {
            this->_super->set(this->_pos, v);
            return *this;
        }
        inline point_reference& operator=(const value_type& v) noexcept(NO_EXCEPT) {
            this->_super->set(this->_pos, v);
            return *this;
        }

        inline point_reference& add(const value_type& v) noexcept(NO_EXCEPT) {
            this->_super->add(this->_pos, v);
            return *this;
        }
        inline point_reference& operator+=(const value_type& v) noexcept(NO_EXCEPT) {
            this->_super->add(this->_pos, v);
            return *this;
        }
    };

    struct range_reference : internal::range_reference<segment_tree> {
        range_reference(segment_tree *const super, const size_type l, const size_type r) noexcept(NO_EXCEPT)
          : internal::range_reference<segment_tree>(super, super->_positivize_index(l), super->_positivize_index(r))
        {
            assert(0 <= this->_begin && this->_begin <= this->_end && this->_end <= this->_super->_n);
        }

        inline value_type fold() noexcept(NO_EXCEPT) {
            if(this->_begin == 0 and this->_end == this->_super->_n) return this->_super->fold();
            return this->_super->fold(this->_begin, this->_end);
        }
    };


    inline auto& add(const size_type p, const value_type& x) noexcept(NO_EXCEPT) {
        assert(0 <= p && p < this->_n);
        this->core::add(p, x);
         return *this;
    }

    inline auto& set(const size_type p, const value_type& x) noexcept(NO_EXCEPT) {
        assert(0 <= p && p < this->_n);
        this->core::set(p, x);
         return *this;
    }

    inline value_type get(const size_type p) const noexcept(NO_EXCEPT) {
        assert(0 <= p && p < this->_n);
        return this->core::fold(p, p+1);
    }

    inline point_reference operator[](const size_type p) noexcept(NO_EXCEPT) { return point_reference(this, p); }
    inline range_reference operator()(const size_type l, const size_type r) noexcept(NO_EXCEPT) { return range_reference(this, l, r); }

    inline value_type fold(const size_type l, const size_type r) const noexcept(NO_EXCEPT) {
        assert(0 <= l && l <= r && r <= this->_n);
        return this->core::fold(l, r);
    }
    inline value_type fold(const size_type r) const noexcept(NO_EXCEPT) {
        assert(0 <= r && r <= this->_n);
        return this->core::fold(0, r);
    }
    inline value_type fold() const noexcept(NO_EXCEPT) {
        return this->core::fold_all();
    }


    template<bool (*f)(value_type)>
    inline size_type max_right(const size_type l) const noexcept(NO_EXCEPT) {
        return this->max_right(l, [](value_type x) { return f(x); });
    }

    template<class F>
    inline size_type max_right(size_type l, F f) const noexcept(NO_EXCEPT) {
        assert(0 <= l && l <= this->_n);
        assert(f(value_type{}));

        if(l == this->_n) return this->_n;

        l += this->_size;
        value_type sm;
        do {
            while((l & 1) == 0) l >>= 1;

            if(!f(sm + this->_data[l])) {
                while(l < this->_size) {
                    l <<= 1;

                    if(f(sm + this->_data[l])) {
                        sm = sm + this->_data[l];
                        ++l;
                    }
                }

                return l - this->_size;
            }

            sm = sm + this->_data[l];
            ++l;
        } while((l & -l) != l);

        return this->_n;
    }


    template<bool (*f)(value_type)>
    inline size_type min_left(const size_type r) const noexcept(NO_EXCEPT) {
        return this->min_left(r, [](value_type x) { return f(x); });
    }

    template<class F>
    inline size_type min_left(size_type r, F f) const noexcept(NO_EXCEPT) {
        assert(0 <= r && r <= this->_n);
        assert(f(value_type{}));

        if (r == 0) return 0;

        r += this->_size;
        value_type sm;

        do {
            --r;
            while(r > 1 && (r & 1)) r >>= 1;

            if(!f(this->_data[r] + sm)) {
                while(r < this->_size) {

                    r = (r << 1 | 1);
                    if(f(this->_data[r] + sm)) {
                        sm = this->_data[r] + sm;
                        --r;
                    }
                }

                return r + 1 - this->_size;
            }

            sm = this->_data[r] + sm;
        } while((r & -r) != r);

        return 0;
    }


    struct iterator : internal::container_iterator_interface<value_type, const segment_tree, iterator> {
        iterator() noexcept = default;

        iterator(const segment_tree *const ref, const size_type p) noexcept(NO_EXCEPT)
          : internal::container_iterator_interface<value_type, const segment_tree, iterator>(ref, p)
        {}

        inline value_type operator*() const noexcept(NO_EXCEPT) { return this->ref()->get(this->pos()); }
    };

    inline iterator begin() const noexcept(NO_EXCEPT) { return iterator(this, 0); }
    inline iterator end() const noexcept(NO_EXCEPT) { return iterator(this, this->_n); }
};


template<actions::internal::operatable_action Action>
struct segment_tree<Action> : segment_tree<typename Action::operand> {
    using segment_tree<typename Action::operand>::segment_tree;
};


} // namespace lib
