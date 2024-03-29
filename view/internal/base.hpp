#pragma once


#include <type_traits>

#include "internal/type_traits.hpp"


namespace lib {

namespace internal {

namespace view_impl {


struct base {};

struct iterator_base {};


template<class T> using is_view = std::is_base_of<base,T>;
template<class T> inline constexpr auto is_view_v = is_view<T>::value;
template<class T> using is_view_t = std::enable_if_t<is_view_v<T>>;

template<class T> using is_view_iterator = std::is_base_of<iterator_base,T>;
template<class T> inline constexpr auto is_view_iterator_v = is_view_iterator<T>::value;
template<class T> using is_view_iterator_t = std::enable_if_t<is_view_iterator_v<T>>;


} // namespace view_impl

} // namespace internal

} // namespace lib
