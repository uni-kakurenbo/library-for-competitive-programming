#pragma once


#include "data_structure/internal/declarations.hpp"

#include "data_structure/range_action/base.hpp"
#include "data_structure/range_action/flags.hpp"

#include "data_structure/monoid/maximum.hpp"
#include "data_structure/monoid/assignment.hpp"


namespace lib {

namespace actions {


template<class T> struct range_set_range_max : base<monoids::assignment<T>>  {
    static constexpr flags tags{ flags::implicit_treap, flags::lazy_segment_tree };

    using operand = monoids::maximum<T>;
    using operation = monoids::assignment<T>;

    static operand map(const operand& x, const operation& y) { return y->value_or(x.val()); }
};


} // namespace actions


template<class T> struct lazy_segment_tree<actions::range_set_range_max<T>> : internal::lazy_segment_tree_lib::core<actions::range_set_range_max<T>> {
    using internal::lazy_segment_tree_lib::core<actions::range_set_range_max<T>>::core;

    inline auto set(const size_t first, const size_t last, const T& val) { return this->apply(first, last, val); }
    inline auto set(const size_t pos, const T& val) { return this->apply(pos, val); }
    inline auto set(const T& val) { return this->apply(val); }

    inline auto max(const size_t first, const size_t last) { return this->fold(first, last); }
    inline auto max() { return this->fold(); }
};


} // namespace lib
