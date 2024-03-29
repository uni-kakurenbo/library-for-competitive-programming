#pragma once


#include <cassert>
#include <type_traits>

#include "internal/dev_env.hpp"


namespace lib {


template<class Engine>
struct random_adaptor {
    using result_type = typename Engine::result_type;
    using signed_result_type = typename std::make_signed_t<result_type>;

  private:
    Engine engine;

  public:
    static constexpr result_type MIN = Engine::min();
    static constexpr result_type MAX = Engine::max();

    static constexpr result_type min() noexcept(NO_EXCEPT) { return MIN; }
    static constexpr result_type max() noexcept(NO_EXCEPT) { return MAX; }

    constexpr random_adaptor(unsigned long seed = 3141592653UL) noexcept(NO_EXCEPT) { this->engine.seed(seed); };

    inline constexpr result_type operator()() noexcept(NO_EXCEPT) {
        return this->engine();
    }

    inline constexpr result_type operator()(result_type max) noexcept(NO_EXCEPT) {
        if(max == 0) return 0;
        return (*this)() % max;
    }
    inline constexpr signed_result_type operator()(const signed_result_type min, const signed_result_type max) noexcept(NO_EXCEPT) {
        assert(min <= max);
        return min + (*this)(max - min);
    };

    template<class T = double>
    inline constexpr T real() noexcept(NO_EXCEPT) {
        const T v = static_cast<T>((this->engine() + 0.5) / (1.0 + this->max()));
        return static_cast<T>((this->operator()() + v) / (1.0 + this->max()));
    }
};


} // namespace lib
