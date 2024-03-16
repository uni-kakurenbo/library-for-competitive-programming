#pragma once


#include <cassert>
#include <vector>
#include <array>

#include <concepts>


#include "internal/dev_env.hpp"

#include "internal/exception.hpp"
#include "internal/types.hpp"

#include "adapter/valarray.hpp"
#include "adapter/vector.hpp"


namespace lib {

namespace internal {

namespace multi_container_impl {


template<class Holder> struct base : Holder {
    using Holder::Holder;

  protected:
    template<std::integral T>
    inline auto _positivize_index(const T _x) const noexcept(NO_EXCEPT) {
        auto x = static_cast<internal::size_t>(_x);
        return x < 0 ? this->size() + x : x;
    }
};


} // namespace multi_contatiner_impl

} // namespace internal


template<class T, unsigned int RANK, template<class...> class Holder = vector, template<class...> class Container = valarray>
struct multi_container : internal::multi_container_impl::base<Holder<multi_container<T, RANK - 1, Holder, Container>>> {
  private:
    using base = internal::multi_container_impl::base<Holder<multi_container<T, RANK - 1, Holder, Container>>>;

  public:
    using base::base;

    template<std::integral Head, class... Tail>
    multi_container(const Head head, Tail&&... tail) noexcept(NO_EXCEPT)
    : base(head, multi_container<T, RANK - 1, Holder, Container>(std::forward<Tail>(tail)...)) {
        assert(head >= 0);
    }

    template<std::integral Head, class... Tail>
    T& operator()(Head _head, Tail&&... tail) noexcept(NO_EXCEPT) {
        static_assert(std::is_integral_v<Head>, "index must be integral");

        const auto index = this->_positivize_index(_head);
        assert(0 <= index && index < std::ranges::ssize(*this));

        return (*this)[index](std::forward<Tail>(tail)...);
    }

    template<std::integral Head, class... Tail>
    const T& operator()(Head _head, Tail&&... tail) const noexcept(NO_EXCEPT) {
        static_assert(std::is_integral_v<Head>, "index must be integral");

        const auto index = this->_positivize_index(_head);
        assert(0 <= index && index < std::ranges::ssize(*this));

        return (*this)[index](std::forward<Tail>(tail)...);
    }
};

template<class T, template<class...> class Holder, template<class...> class Container>
struct multi_container<T, 1, Holder, Container> : internal::multi_container_impl::base<Container<T>> {
    using internal::multi_container_impl::base<Container<T>>::base;

    template<class... Args>
    multi_container(const Args&... args) noexcept(NO_EXCEPT) : internal::multi_container_impl::base<Container<T>>(args...)
    {}

    template<class Index>
    T& operator()(Index&& _index) noexcept(NO_EXCEPT) {
        const auto index = this->_positivize_index(std::forward<T>(_index));
        assert(0 <= index && index < std::ranges::ssize(*this));

        return (*this)[index];
    }

    template<class Index>
    const T& operator()(Index&& _index) const noexcept(NO_EXCEPT) {
        const auto index = this->_positivize_index(std::forward<T>(_index));
        assert(0 <= index && index < std::ranges::ssize(*this));

        return (*this)[index];
    }
};


template<class T, template<class...> class Holder, template<class...> class Container>
struct multi_container<T, 0, Holder, Container> {
    static_assert(internal::EXCEPTION<T>, "invalid rank: 0, should be 1 or more");
};


} // namespace lib
