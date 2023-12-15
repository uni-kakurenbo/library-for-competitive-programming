#pragma once


#include <cassert>
#include <cstdint>
#include <iostream>


#include "snippet/aliases.hpp"

#include "internal/dev_env.hpp"
#include "internal/types.hpp"

#include "numeric/modular/modint_interface.hpp"
// #include "numeric/internal/primality_test.hpp"
// #include "numeric/internal/primitive_root.hpp"

#include "template/debug.hpp"


namespace lib {


template<internal::modular_context Context>
struct modint {
    using value_type = typename Context::value_type;

  private:
    using mint = modint;
    using context = Context;

    value_type _val = 0;


    static constexpr mint _raw(const value_type v) noexcept(NO_EXCEPT) {
        mint res;
        res._val = v;
        return res;
    }


    static constexpr mint _one() noexcept(NO_EXCEPT)
        requires internal::has_static_one<typename mint::context::reduction>
    {
        return mint::context::get().one;
    }

    static constexpr mint _one() noexcept(NO_EXCEPT)
        requires (!internal::has_static_one<typename mint::context::reduction>)
    {
        return 1;
    }

  public:
    static constexpr int digits = mint::context::reduction::digits;
    static inline constexpr value_type max() noexcept { return mint::context::reduction::max(); }

    static constexpr mint zero = mint::_raw(0);
    static inline mint one = mint::_one();


    static constexpr void set_mod(const value_type mod) noexcept(NO_EXCEPT)
        requires requires(value_type mod) { mint::context::set_mod(mod); }
    {
        mint::context::set_mod(mod);
        mint::one = mint::_one();
    }


    static inline constexpr value_type mod() noexcept(NO_EXCEPT) { return mint::context::get().mod(); }

    static inline constexpr mint raw(const value_type v) noexcept(NO_EXCEPT)
    {
        mint res;
        res._val = mint::context::get().convert_raw(v);
        return res;
    };


    constexpr modint() noexcept = default;

    template<std::integral T>
    constexpr modint(const T v) noexcept(NO_EXCEPT) : _val(mint::context::get().convert(v)) {}


    inline constexpr value_type val() const noexcept(NO_EXCEPT) {
        return mint::context::get().revert(this->_val);
    }

    inline constexpr explicit operator value_type() const noexcept(NO_EXCEPT) { return this->val(); }


    inline constexpr mint& operator+=(const mint& rhs) noexcept(NO_EXCEPT) {
        this->_val = mint::context::get().add(this->_val, rhs._val);
        return *this;
    }

    inline constexpr mint& operator-=(const mint& rhs) noexcept(NO_EXCEPT) {
        this->_val = mint::context::get().subtract(this->_val, rhs._val);
        return *this;
    }


    inline constexpr mint& operator*=(const mint& rhs) noexcept(NO_EXCEPT) {
        this->_val = mint::context::get().multiply(this->_val, rhs._val);
        return *this;
    }

    inline constexpr auto& operator/=(const mint& rhs) noexcept(NO_EXCEPT) { return *this *= rhs.inv(); }


    constexpr mint pow(const i64 n) const noexcept(NO_EXCEPT) {
        return mint::_raw(mint::context::get().pow(this->_val, n));
    }


    constexpr mint inv() const noexcept(NO_EXCEPT) {
        using signed_value_type = std::make_signed_t<value_type>;

        signed_value_type x = this->val(), y = mint::mod(), u = 1, v = 0;

        while(y > 0) {
            signed_value_type t = x / y;
            std::swap(x -= t * y, y);
            std::swap(u -= t * v, v);
        }

        assert(x == 1);

        if(u < 0) u += v / x;
        return mint::raw(u);
    }

    friend inline constexpr bool operator==(const mint& lhs, const mint& rhs) noexcept(NO_EXCEPT) {
        return mint::context::get().equal(lhs._val, rhs._val);
    }

    friend inline constexpr bool operator!=(const mint& lhs, const mint& rhs) noexcept(NO_EXCEPT) { return !(lhs == rhs); }


    inline constexpr mint& operator++() noexcept(NO_EXCEPT) { return *this += mint::one; }
    inline constexpr mint& operator--() noexcept(NO_EXCEPT) { return *this -= mint::one; }

    inline constexpr mint operator++(int) noexcept(NO_EXCEPT) { const mint res = *this; return ++*this, res; }
    inline constexpr mint operator--(int) noexcept(NO_EXCEPT) { const mint res = *this; return --*this, res; }

    inline constexpr auto operator+() const noexcept(NO_EXCEPT) { return *this; }
    inline constexpr auto operator-() const noexcept(NO_EXCEPT) { return mint::zero - *this; }

    friend inline constexpr mint operator+(mint lhs, const mint& rhs) noexcept(NO_EXCEPT) { return lhs += rhs; }
    friend inline constexpr mint operator-(mint lhs, const mint& rhs) noexcept(NO_EXCEPT) { return lhs -= rhs; }
    friend inline constexpr mint operator*(mint lhs, const mint& rhs) noexcept(NO_EXCEPT) { return lhs *= rhs; }
    friend inline constexpr mint operator/(mint lhs, const mint& rhs) noexcept(NO_EXCEPT) { return lhs /= rhs; }
};


} // namespace lib