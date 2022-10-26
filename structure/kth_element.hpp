#include <vector>
#include <functional>
#include <queue>

template<class T, class Container = std::vector<T>, class CompA = std::less<T>, class CompB = std::greater<T>>
struct KthElement {
  private:
    size_t K;
    std::priority_queue<T,Container,CompA> small;
    std::priority_queue<T,Container,CompB> large;

  public:
    KthElement(size_t K) : K(K) {}

    inline T get() const { return small.top(); }
    inline bool has() const { return small.size() == K; }

    inline void push(T v) {
        if(small.size() < K) {
            small.push(v);
            return;
        }
        T kth = small.top();
        if(v < kth) {
            small.pop(); small.push(v);
            large.push(kth);
        }
        else {
            large.push(v);
        }
    }

    inline void pop() {
        small.pop();
        if(large.empty()) return;
        T v = large.top(); large.pop();
        small.push(v);
    }
};
