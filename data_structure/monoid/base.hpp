#pragma once

#include "data_structure/semigroup/base.hpp"

namespace lib {

namespace monoids {

namespace internal {

struct base : lib::semigroups::internal::base {};

} // namespace internal


template<class T> struct base : internal::base {
    using value_type = T;

  protected:
    value_type _value;

  public:
    base(const value_type& value = {}) : _value(value) {};

    inline operator value_type() const { return this->_value; }
    inline value_type val() const { return this->_value; };

    inline const value_type* operator->() const { return &this->_value; };
    inline value_type* operator->() { return &this->_value; };

    inline bool operator==(const base& other) const { return this->_value == other._value; };
    inline bool operator!=(const base& other) const { return this->_value != other._value; };
};


} // namespace internal

} // namespace lib
