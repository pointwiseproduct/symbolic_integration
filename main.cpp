//-------- include start --------//
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <sstream>
#include <regex>
#include <functional>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <cstdlib>
#include "mpirxx.h"

namespace symbolic_alg{

class poly;
poly operator +(const poly&, const poly&);
poly operator -(const poly&, const poly&);
poly operator *(const poly&, const poly&);

namespace aux{
    template<class T>
    struct rebind_1st;

    template<class P, template<class> class T>
    struct rebind_1st<T<P>>{
        template<class U>
        struct rebind{
            using type = T<U>;
        };
    };

    template<
        class Key,
        class Mapped,
        class Compare = std::less<Key>,
        class Alloc = std::allocator<std::pair<Key, Mapped>>
    > class vector_map{
    public:
        using key_type = Key;
        using mapped_type = Mapped;
        using key_compare = Compare;
        using value_type = std::pair<key_type, mapped_type>;
        using container_type = std::vector<std::pair<key_type, mapped_type>>;
        using allocator_type = Alloc;
        using reference = value_type&;
        using const_reference = const value_type&;
        using pointer = typename std::allocator_traits<allocator_type>::pointer;
        using const_pointer = typename std::allocator_traits<allocator_type>::const_pointer;
        using iterator = typename container_type::iterator;
        using const_iterator = typename container_type::const_iterator;
        using reverse_iterator = typename container_type::reverse_iterator;
        using const_reverse_iterator = typename container_type::const_reverse_iterator;
        using difference_type = typename container_type::difference_type;
        using size_type = typename container_type::size_type;

        class value_compare{
            friend class vector_map;
        public:
            Compare comp;
            value_compare(Compare comp) : comp(comp){}
            using result_type = bool;
            using first_argument_type = value_type;
            using second_argument_type = value_type;
            bool operator ()(const value_type& x, const value_type& y) const{
                return comp(x.first, y.first);
            }
        };

        explicit vector_map(
            const key_compare &comp = key_compare(),
            const allocator_type &alloc = allocator_type()
        ) : comp(comp), vec(alloc){}

        explicit vector_map(const allocator_type &alloc) : comp(), vec(alloc){}

        template<class InputIterator>
        vector_map(
            InputIterator first, InputIterator last,
            const key_compare &comp = key_compare(),
            const allocator_type &alloc = allocator_type()
        ) : comp(comp), vec(first, last, alloc)
        { std::sort(vec.begin(), vec.end(), comp); }

        vector_map(const vector_map &other)
            : comp(other.comp), vec(other.vec)
        {}

        vector_map(vector_map &&other)
            : comp(std::move(other.comp)), vec(std::move(other.vec))
        {}

        vector_map(vector_map &&other, const allocator_type &alloc)
            : comp(std::move(otehr.comp)), vec(std::move(otehr.vec), alloc)
        {}

        vector_map(
            std::initializer_list<value_type> list,
            const key_compare &comp = key_compare(),
            const allocator_type &alloc = allocator_type()
        ) : comp(comp), vec(list, alloc)
        { std::sort(vec.begin(), vec.end(), [&](const value_type &a, const value_type &b){ return comp(a.first, b.first); }); }

        vector_map &operator =(const vector_map &other){
            comp = other.comp;
            vec = other.vec;
            return *this;
        }

        vector_map &operator =(vector_map &&other){
            comp = std::move(other.comp);
            vec = std::move(other.vec);
            return *this;
        }

        vector_map &operator =(std::initializer_list<value_type> list){
            vec = list;
            return *this;
        }

        iterator begin(){
            return vec.begin();
        }

        iterator end(){
            return vec.end();
        }

        const_iterator begin() const{
            return vec.begin();
        }

        const_iterator end() const{
            return vec.end();
        }

        iterator rbegin(){
            return vec.rbegin();
        }

        iterator rend(){
            return vec.rend();
        }

        const_iterator rbegin() const{
            return vec.rbegin();
        }

        const_iterator rend() const{
            return vec.rend();
        }

        const_iterator cbegin() const{
            return vec.cbegin();
        }

        const_iterator cend() const{
            return vec.cend();
        }

        const_iterator crbegin() const{
            return vec.crbegin();
        }

        const_iterator crend() const{
            return vec.crend();
        }

        bool empty() const{
            return vec.empty();
        }

        size_type size() const{
            return vec.size();
        }

        bool max_size() const{
            return vec.max_size();
        }

        value_type &operator [](const key_type &key){
            std::pair<key_type, mapped_type> v;
            v.first = key;
            auto iter = std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == key){
                return *iter;
            }else{
                return vec.end();
            }
        }

        const value_type &operator [](const key_type &key) const{
            std::pair<key_type, mapped_type> v;
            v.first = key;
            auto iter = std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == key){
                return *iter;
            }else{
                return vec.end();
            }
        }

        value_type &at(const key_type &key){
            std::pair<key_type, mapped_type> v;
            v.first = key;
            auto iter = std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == key){
                return *iter;
            }else{
                throw std::out_of_range("vector_map: throw out_of_range;");
            }
        }

        const value_type &at(const key_type &key) const{
            std::pair<key_type, mapped_type> v;
            v.first = key;
            auto iter = std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == key){
                return *iter;
            }else{
                throw std::out_of_range("vector_map: throw out_of_range;");
            }
        }

        std::pair<iterator, bool> insert(const value_type &value){
            iterator iter = std::lower_bound(
                vec.begin(), vec.end(), value,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == value.first){
                return std::make_pair(iter, false);
            }else{
                iter = vec.insert(iter, value);
                return std::make_pair(iter, true);
            }
        }

        std::pair<iterator, bool> insert(value_type &&value){
            iterator iter = std::lower_bound(
                vec.begin(), vec.end(), value,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == value.first){
                return std::make_pair(iter, false);
            }else{
                iter = vec.insert(iter, value);
                return std::make_pair(iter, true);
            }
        }

        iterator insert(const_iterator position, const value_type &value){
            return vec.insert(position, value);
        }

        iterator insert(const_iterator position, value_type &&value){
            return vec.insert(position, value);
        }

        template<class InputIter>
        void insert(InputIter first, InputIter last){
            vec.insert(vec.end(), first, last);
            std::sort(
                vec.begin(), vec.end(),
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
        }

        void insert(std::initializer_list<value_type> list){
            vec.insert(list);
            std::sort(
                vec.begin(), vec.end(),
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
        }

        iterator erase(const_iterator position){
            return vec.erase(position);
        }

        iterator erase(const_iterator first, const_iterator last){
            return vec.erase(first, last);
        }

        size_type erase(const key_type &key){
            std::pair<key_type, mapped_type> v;
            v.first = key;
            iterator iter = std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == key){
                vec.erase(iter);
                return 1;
            }else{
                return 0;
            }
        }

        void swap(vector_map &other){
            std::swap(comp, other.comp);
            vec.swap(other.vec);
        }

        void clear() noexcept{
            vec.clear();
        }

        std::pair<iterator, bool> emplace(key_type &&key, mapped_type &&mapped){
            std::pair<key_type, mapped_type> v = { key, mapped };
            iterator iter = std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == value.first){
                return std::make_pair(iter, false);
            }else{
                iter = vec.insert(iter, v);
                return std::make_pair(iter, true);
            }
        }

        std::pair<iterator, bool> emplace(const_iterator hint, key_type &&key, mapped_type &&mapped){
            std::pair<key_type, mapped_type> v = { key, mapped };
            if(hint == end()){
                --hint;
            }
            if(comp(*hint, v)){
                return insert(v);
            }else if(!comp(v, *hint)){
                return std::make_pair(hint, false);
            }else{
                iterator iter = vec.insert(hint, v);
                return std::make_pair(iter, true);
            }
        }

        key_compare key_comp() const{
            return key_compare();
        }

        value_compare value_comp() const{
            return value_compare(comp);
        }

        iterator find(const key_type &key){
            std::pair<key_type, mapped_type> v;
            v.first = key;
            iterator iter = std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == key){
                return iter;
            }else{
                return vec.end();
            }
        }

        const_iterator find(const key_type &key) const{
            std::pair<key_type, mapped_type> v;
            v.first = key;
            const_iterator iter = std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == key){
                return iter;
            }else{
                return vec.end();
            }
        }

        size_type count(const key_type &key) const{
            if(find(key) != end()){
                return 1;
            }else{
                return 0;
            }
        }

        iterator lower_bound(const key_type &key){
            std::pair<key_type, mapped_type> v;
            v.first = key;
            return std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
        }

        const_iterator lower_bound(const key_type &key) const{
            std::pair<key_type, mapped_type> v;
            v.first = key;
            return std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
        }

        iterator upper_bound(const key_type &key){
            std::pair<key_type, mapped_type> v;
            v.first = key;
            return std::upper_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
        }

        const_iterator upper_bound(const key_type &key) const{
            std::pair<key_type, mapped_type> v;
            v.first = key;
            return std::upper_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
        }

        std::pair<iterator, iterator> equal_range(const key_type &key){
            std::pair<key_type, mapped_type> v;
            v.first = key;
            iterator iter = std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == key){
                std::make_pair(iter, iter);
            }
        }

        std::pair<const_iterator, const_iterator> equal_range(const key_type &key) const{
            std::pair<key_type, mapped_type> v;
            v.first = key;
            const_iterator iter = std::lower_bound(
                vec.begin(), vec.end(), v,
                [&](const value_type &a, const value_type &b){
                    return comp(a.first, b.first);
                }
            );
            if(iter != vec.end() && iter->first == key){
                std::make_pair(iter, iter);
            }
        }

        static int lexicographical_key_type_compare(const vector_map &lhs, const vector_map &rhs){
            auto iter = lhs.begin();
            auto jter = rhs.begin();
            while(iter != lhs.end() && jter != rhs.end()){
                if(lhs.comp(iter->first, jter->first)){
                    return -1;
                }else if(!lhs.comp(iter->first, jter->first) && !lhs.comp(jter->first, iter->first)){
                    ++iter, ++jter;
                    continue;
                }else{
                    return +1;
                }
            }
            if(iter == lhs.end() && jter != rhs.end()){
                return -1;
            }else if(iter != lhs.end() && jter == rhs.end()){
                return +1;
            }else{
                return lexicographical_mapped_type_compare(lhs, rhs);
            }
        }

        static int lexicographical_mapped_type_compare(const vector_map &lhs, const vector_map &rhs){
            using comparetor = typename rebind_1st<key_compare>::rebind<mapped_type>::type;
            auto iter = lhs.begin();
            auto jter = rhs.begin();
            while(iter != lhs.end() && jter != rhs.end()){
                if(comparetor()(iter->second, jter->second)){
                    return -1;
                }else if(!comparetor()(iter->second, jter->second) && !comparetor()(jter->second, iter->second)){
                    ++iter, ++jter;
                    continue;
                }else{
                    return +1;
                }
            }
            if(iter == lhs.end() && jter != rhs.end()){
                return -1;
            }else if(iter != lhs.end() && jter == rhs.end()){
                return +1;
            }else{
                return 0;
            }
        }

        bool operator ==(const vector_map &rhs) const{
            return lexicographical_key_type_compare(*this, rhs) == 0;
        }

        bool operator !=(const vector_map &rhs) const{
            return !(*this == rhs);
        }

        bool operator <(const vector_map &rhs) const{
            return lexicographical_key_type_compare(*this, rhs) < 0;
        }

        bool operator <=(const vector_map &rhs) const{
            return lexicographical_key_type_compare(*this, rhs) <= 0;
        }

        bool operator >(const vector_map &rhs) const{
            return lexicographical_key_type_compare(*this, rhs) > 0;
        }

        bool operator >=(const vector_map &rhs) const{
            return lexicographical_key_type_compare(*this, rhs) >= 0;
        }

    private:
        key_compare comp;
        std::vector<std::pair<key_type, mapped_type>> vec;
    };

    template<class T>
    static T pow(T a, int n){
        if(n == 0){
            return T(1);
        }
        bool positive = n > 0;
        if(!positive){
            n = -n;
        }
        T r = T(1);
        for(; n; n >>= 1, a *= a){
            if(n & 1){
                r *= a;
            }
        }
        if(positive){
            return std::move(r);
        }else{
            return 1 / r;
        }
    }

    template<class T>
    static T abs(const T &a){
        return a < 0 ? -a : a;
    }

    struct parsing_semantic_data;
    poly parse(const std::string &str);
}

class poly{
    friend aux::parsing_semantic_data;

public:
    using degree_type = int;

    // 係数の型
    using coefficient_type = mpq_class;

    // 項の型
    class term_type{
    public:
        term_type() = default;
        term_type(degree_type deg, coefficient_type coe)
            : deg(deg), coe(coe)
        {}

        ~term_type() = default;

        degree_type deg;
        coefficient_type coe;
    };

private:
    class term_compare{
    public:
        bool operator()(const term_type &lhs, const term_type &rhs) const{
            return lhs.deg < rhs.deg;
        }
    };

    using data_type = std::vector<term_type>;
    data_type data;

public:
    using iterator = typename data_type::iterator;
    using const_iterator = typename data_type::const_iterator;
    using reverse_iterator = typename data_type::reverse_iterator;
    using const_reverse_iterator = typename data_type::const_reverse_iterator;

    iterator begin(){
        return data.begin();
    }

    const_iterator begin() const{
        return data.begin();
    }

    reverse_iterator rbegin(){
        return data.rbegin();
    }

    const_reverse_iterator rbegin() const{
        return data.rbegin();
    }

    iterator end(){
        return data.end();
    }

    const_iterator end() const{
        return data.end();
    }

    reverse_iterator rend(){
        return data.rend();
    }

    const_reverse_iterator rend() const{
        return data.rend();
    }

    degree_type deg() const{
        if(data.empty()){
            return degree_type(0);
        }
        const_reverse_iterator iter = rbegin();
        return iter->deg;
    }

    coefficient_type lc() const{
        if(data.empty()){
            return coefficient_type(0);
        }
        const_reverse_iterator iter = rbegin();
        return iter->coe;
    }

    poly() = default;
    poly(const poly&) = default;
    poly(poly &&other) : data(std::move(other.data)){}
    poly(const degree_type &d, const coefficient_type &c) : data(){
        if(c != 0){
            data.push_back(term_type(d, c));
        }
    }

    poly(const term_type &term) : data(){
        data.push_back(term);
    }

    poly(const coefficient_type &coe) : data(){
        if(coe != 0){
            data.push_back(term_type(degree_type(0), coe));
        }
    }

    poly(int n) : data(){
        if(n != 0){
            data.push_back(term_type(degree_type(0), coefficient_type(n)));
        }
    }

    ~poly() = default;

    coefficient_type operator ()(const degree_type &i) const{
        term_type t = { i, coefficient_type(0) };
        auto iter = std::lower_bound(data.begin(), data.end(), t, term_compare());
        if(iter->deg == i){
            return iter->coe;
        }else{
            return 0;
        }
    }

    bool is_zero() const{
        return data.empty();
    }

    bool is_not_zero() const{
        return !data.empty();
    }

    void negative(){
        for(auto &&n : data){
            n.coe = -n.coe;
        }
    }

    poly operator -() const{
        poly p = *this;
        p.negative();
        return p;
    }

    void add(const poly &other){
        if(data.empty()){
            *this = other;
            return;
        }
        for(auto &n : other.data){
            auto iter = std::lower_bound(data.begin(), data.end(), n, term_compare());
            if(iter != data.end() && iter->deg == n.deg){
                iter->coe += n.coe;
                if(iter->coe == coefficient_type(0)){
                    data.erase(iter);
                }
            }else{
                data.insert(iter, n);
            }
        }
    }

    void sub(const poly &other){
        if(data.empty()){
            *this = other;
            negative();
            return;
        }
        for(auto &n : other.data){
            auto iter = std::lower_bound(data.begin(), data.end(), n, term_compare());
            if(iter != data.end() && iter->deg == n.deg){
                iter->coe -= n.coe;
                if(iter->coe == coefficient_type(0)){
                    data.erase(iter);
                }
            }else{
                term_type t = n;
                t.coe = -n.coe;
                data.insert(iter, t);
            }
        }
    }

    poly &operator =(const poly &other){
        data = other.data;
        return *this;
    }

    poly &operator =(poly &&other){
        data = std::move(other.data);
        return *this;
    }

    poly &operator +=(const poly &other){
        add(other);
        return *this;
    }

    poly &operator +=(const term_type &t){
        if(data.empty()){
            data.push_back(t);
            return *this;
        }
        auto iter = std::lower_bound(data.begin(), data.end(), t, term_compare());
        if(iter != data.end() && iter->deg == t.deg){
            iter->coe += t.coe;
            if(iter->coe == coefficient_type(0)){
                data.erase(iter);
            }
        }else{
            data.insert(iter, t);
        }
        return *this;
    }

    poly &operator -=(const term_type &t){
        if(data.empty()){
            term_type u = t;
            u.coe = -u.coe;
            data.push_back(u);
            return *this;
        }
        auto iter = std::lower_bound(data.begin(), data.end(), t, term_compare());
        if(iter != data.end() && iter->deg == t.deg){
            iter->coe -= t.coe;
            if(iter->coe == coefficient_type(0)){
                data.erase(iter);
            }
        }else{
            term_type u = t;
            u.coe = -u.coe;
            data.insert(iter, u);
        }
        return *this;
    }

    poly &operator -=(const poly &other){
        sub(other);
        return *this;
    }

    void primitive_mul(const poly &lhs, const poly &rhs){
        *this = coefficient_type(0);
        for(const auto &i : lhs){
            for(const auto &j : rhs){
                *this += term_type(i.deg + j.deg, i.coe * j.coe);
            }
        }
    }

    static poly mul(const poly &lhs, const poly &rhs){
        poly r;
        r.primitive_mul(lhs, rhs);
        return r;
    }

    static std::tuple<poly, poly> euclidean_div(const poly &A, const poly &B){
        poly Q = 0, R = A;
        degree_type delta;
        while(R.is_not_zero() && (delta = R.deg() - B.deg()) >= 0){
            poly T(delta, R.lc() / B.lc());
            Q += T;
            R -= B * T;
        }
        return std::make_tuple(Q, R);
    }

    static std::tuple<poly, poly> psude_div(const poly &A, const poly &B){
        poly b = B.lc();
        degree_type N = A.deg() - B.deg() + 1;
        poly Q = 0, R = A;
        degree_type delta;
        while(R.is_not_zero() && (delta = R.deg() - B.deg()) >= 0){
            poly T(delta, R.lc());
            --N;
            Q = b * Q + T;
            R = b * R - T * B;
        }
        poly bN = pow(b, N);
        return std::make_tuple(bN * Q, bN * R);
    }

    static std::tuple<poly, poly> quo_rem(const poly &A, const poly &B){
        return euclidean_div(A, B);
    }

    static poly quo(const poly &A, const poly &B){
        return std::get<0>(euclidean_div(A, B));
    }

    static poly div(const poly &A, const poly &B){
        auto p = euclidean_div(A, B);
        if(std::get<1>(p) != 0){
            throw;
        }
        return std::get<0>(p);
    }

    poly &operator *=(const poly &other){
        *this = mul(*this, other);
        return *this;
    }

    static int lexicographical_compare(const poly &lhs, const poly &rhs){
        auto iter = lhs.data.rbegin();
        auto jter = rhs.data.rbegin();
        while(iter != lhs.data.rend() && jter != rhs.data.rend()){
            if(iter->deg == jter->deg){
                if(iter->coe == iter->coe){
                    ++iter, ++jter;
                    continue;
                }else if(iter->coe > jter->coe){
                    return -1;
                }else{
                    return +1;
                }
            }else if(iter->deg > jter->deg){
                if(iter->coe < 0){
                    return -1;
                }else{
                    return +1;
                }
            }else{
                if(jter->coe < 0){
                    return +1;
                }else{
                    return -1;
                }
            }
        }
        if(iter == lhs.data.rend() && jter == rhs.data.rend()){
            return 0;
        }else if(iter == lhs.data.rend()){
            return -1;
        }else{
            return +1;
        }
    }

    bool operator ==(const poly &other) const{
        return lexicographical_compare(*this, other) == 0;
    }

    bool operator !=(const poly &other) const{
        return !(*this == other);
    }

    class content_exception : public std::runtime_error{
    public:
        content_exception() : std::runtime_error("coefficient multiplies is not integer."){}
        content_exception(const content_exception&) = default;
    };

    coefficient_type content() const{
        if(is_zero()){
            return 0;
        }
        bool all_int = true;
        for(const term_type &t : data){
            if(t.coe.get_den() != 1){
                all_int = false;
                break;
            }
        }
        if(all_int){
            int b = coe_to_int(data[0].coe);
            for(std::size_t i = 1; i < data.size(); ++i){
                b = gcd(b, coe_to_int(data[i].coe));
            }
            return b;
        }else{
            mpz_class n = abs(data[0].coe.get_num()), d = abs(data[0].coe.get_den());
            for(std::size_t i = 1; i < data.size(); ++i){
                if(!(data[i].coe.get_num() % n == 0 && data[i].coe.get_den() % d == 0)){
                    bool divid = false;
                    if(n % data[i].coe.get_num() == 0){
                        n = abs(data[i].coe.get_num());
                        divid = true;
                    }
                    if(d % data[i].coe.get_den() == 0){
                        d = abs(data[i].coe.get_den());
                        divid = true;
                    }
                    if(!divid){
                        throw content_exception();
                    }
                }
            }
            return coefficient_type(n) / coefficient_type(d);
        }
    }

    poly pp() const{
        if(is_zero()){
            return 0;
        }
        coefficient_type c = content();
        poly r = *this;
        for(auto &t : r.data){
            t.coe /= c;
        }
        return r;
    }

    poly diff() const{
        return diff(1);
    }

    poly diff(std::size_t n) const{
        poly p = *this;
        for(term_type &t : p){
            t.coe *= t.deg;
            --t.deg;
        }
        for(std::size_t i = 0; i < p.data.size(); ++i){
            auto iter = p.data.begin() + i;
            if(iter->coe == 0){
                p.data.erase(iter);
                --i;
            }
        }
        return p;
    }

    static poly euclidean(poly a, poly b){
        if(a.deg() < b.deg()){
            std::swap(a, b);
        }
        while(b.is_not_zero()){
            auto p = euclidean_div(a, b);
            a = b;
            b = std::get<1>(p);
        }
        return a;
    }

    static poly gcd(poly a, poly b){
        return euclidean(a, b).pp();
    }

    static std::tuple<poly, poly, poly> extended_euclidean(poly a, poly b){
        poly a1 = 1, a2 = 0, b1 = 0, b2 = 1;
        while(b.is_not_zero()){
            auto qr = euclidean_div(a, b);
            a = b;
            b = std::get<1>(qr);
            poly r1 = a1 - std::get<0>(qr) * b1, r2 = a2 - std::get<0>(qr) * b2;
            a1 = b1;
            a2 = b2;
            b1 = r1;
            b2 = r2;
        }
        return std::make_tuple(a1, a2, a);
    }

    static std::tuple<poly, poly> half_extended_euclidean(poly a, poly b){
        poly a1 = 1, b1 = 0;
        while(b.is_not_zero()){
            auto qr = quo_rem(a, b);
            a = b;
            b = std::get<1>(qr);
            poly r1 = a1 - std::get<0>(qr) * b1;
            a1 = b1;
            b1 = r1;
        }
        return std::make_tuple(a1, a);
    }

    static std::tuple<poly, poly, poly> half_full_extended_euclidean(const poly &a, const poly &b){
        auto sg = half_extended_euclidean(a, b);
        auto tr = quo_rem(std::get<1>(sg) - std::get<0>(sg) * a, b);
        return std::make_tuple(std::get<0>(sg), std::get<0>(tr), std::get<1>(sg));
    }

    class extended_euclidean_exception : public std::runtime_error{
    public:
        extended_euclidean_exception(const char *msg) : std::runtime_error(msg){}
        extended_euclidean_exception(const extended_euclidean_exception &other) : std::runtime_error(other){}
    };

    static std::tuple<poly, poly> extended_euclidean(const poly &a, const poly &b, const poly &c){
        poly s, t, g, q, r;
        std::tie(s, t, g) = half_full_extended_euclidean(a, b);
        std::tie(q, r) = quo_rem(c, g);
        if(r.is_not_zero()){
            throw extended_euclidean_exception("c is not in the ideal generated by a and b.");
        }
        s = q * s;
        t = q * t;
        if(s.is_not_zero() && s.deg() >= b.deg()){
            std::tie(q, r) = quo_rem(s, b);
            s = r;
            t = t + q * a;
        }
        return std::make_tuple(s, t);
    }

    static poly half_extended_euclidean(const poly &a, const poly &b, const poly &c){
        poly s, g, q, r;
        std::tie(s, g) = half_extended_euclidean(a, b);
        std::tie(q, r) = quo_rem(c, g);
        if(r.is_not_zero()){
            throw extended_euclidean_exception("c is not in the ideal generated by a and b.");
        }
        s = q * s;
        if(s.is_not_zero() && s.deg() >= b.deg()){
            std::tie(q, r) = quo_rem(s, b);
            s = r;
        }
        return s;
    }

    static std::tuple<poly, poly> half_full_extended_euclidean(const poly &a, const poly &b, const poly &c){
        poly s = half_extended_euclidean(a, b, c), t, r;
        std::tie(t, r) = quo_rem(c - s * a, b);
        return std::make_tuple(s, t);
    }

    template<class... T>
    static std::vector<poly> partial_fraction(const poly &a, const poly &d1, const T&... d){
        std::vector<poly> v({ d1, d... });
        return partial_fraction(a, std::move(v));
    }

    static std::vector<poly> partial_fraction(const poly &a, std::vector<poly> d){
        std::vector<poly> result;
        std::size_t n = d.size();
        poly a0, r, d_prod = 1, d_prod2 = 1;
        for(std::size_t i = 1; i < n; ++i){
            d_prod2 *= d[i];
        }
        d_prod = d_prod2 * d[0];
        std::tie(a0, r) = quo_rem(a, d_prod);
        if(n == 1){
            result = { a0, r };
            return result;
        }
        poly a1, t;
        std::tie(a1, t) = extended_euclidean(d_prod2, d[0], r);
        d.erase(d.begin());
        std::vector<poly> rec_result = partial_fraction(t, std::move(d));
        rec_result[0] += a0;
        rec_result.insert(rec_result.begin() + 1, a1);
        return std::move(rec_result);
    }

    static std::vector<poly> partial_fraction(const poly &a, const std::vector<poly> &d, const std::vector<degree_type> &e){
        std::vector<poly> de(d.size());
        for(std::size_t i = 0; i < de.size(); ++i){
            de[i] = pow(d[i], e[i]);
        }
        std::vector<poly> an = partial_fraction(a, de);
        poly a0 = std::move(an[0]);
        an.erase(an.begin());
        std::vector<poly> result;
        for(std::size_t i = 0; i < an.size(); ++i){
            for(degree_type j = e[i]; j > 0; --j){
                poly q, aij;
                std::tie(q, aij) = quo_rem(an[i], d[i]);
                an[i] = q;
                result.push_back(std::move(aij));
            }
            std::reverse(result.end() - e[i], result.end());
            an[0] = an[0] + an[i];
        }
        result.insert(result.begin(), a0);
        return std::move(result);
    }

    static std::tuple<poly, std::vector<poly>> sub_resultant(const poly &A, const poly &B){
        std::vector<poly> R;
        R.push_back(A);
        R.push_back(B);
        std::size_t i = 1;
        std::vector<coefficient_type> γ, β, r;
        std::vector<degree_type> δ;
        γ.push_back(0); γ.push_back(-1);
        δ.push_back(0); δ.push_back(A.deg() - B.deg());
        β.push_back(0); β.push_back(pow(-1, δ[1] + 1));
        r.push_back(0);
        while(R[i].is_not_zero()){
            r.push_back(R[i].lc());
            poly Quo, Rem;
            std::tie(Quo, Rem) = psude_div(R[i - 1], R[i]);
            R.push_back(std::get<0>(quo_rem(Rem, poly(β[i]))));
            ++i;
            γ.push_back(aux::pow(coefficient_type(-γ[i - 1]), δ[i - 1]) * aux::pow(γ[i - 1], 1 - δ[i - 1]));
            δ.push_back(R[i - 1].deg() - R[i].deg());
            β.push_back(-r[i - 1] * pow(γ[i], δ[i]));
        }
        std::size_t k = i - 1;
        R[i] = poly(0);
        if(R[k].deg() > 0){
            return std::make_tuple(poly(0), std::move(R));
        }
        if(R[k - 1].deg() == 1){
            return std::make_tuple(R[k], std::move(R));
        }
        coefficient_type s = 1, c = 1;
        for(std::size_t j = 1; j < k - 1; ++j){
            if(odd(R[j - 1].deg()) && odd(R[j])){
                s = -s;
            }
            c = c * aux::pow(coefficient_type(β[j] / pow(r[j], 1 + δ[j])), R[j].deg()) * aux::pow(r[j], R[j - 1].deg() - R[j + 1].deg());
        }
        return std::make_tuple(poly(s) * poly(c) * pow(R[k], R[k - 1].deg()), std::move(R));
    }

    static std::vector<poly> squarefree_musser(const poly &A){
        coefficient_type c = A.content();
        poly S = std::get<0>(quo_rem(A, poly(c)));
        poly S_minus = gcd(S, S.diff());
        poly S_star = std::get<0>(quo_rem(S, S_minus));
        std::size_t k = 0;
        std::vector<poly> Ak;
        while(S_minus.deg() > 0){
            poly Y = gcd(S_star, S_minus);
            Ak.push_back(std::get<0>(quo_rem(S_star, Y)));
            S_star = Y;
            S_minus = std::get<0>(quo_rem(S_minus, Y));
            ++k;
        }
        Ak.push_back(S_star);
        Ak[0] *= poly(c) * S_minus;
        return std::move(Ak);
    }

    static std::vector<poly> squarefree_yun(const poly &A){
        coefficient_type c = A.content();
        poly S = std::get<0>(quo_rem(A, poly(c)));
        poly S_prime = S.diff();
        poly S_minus = gcd(S, S_prime);
        poly S_star = std::get<0>(quo_rem(S, S_minus));
        poly Y = std::get<0>(quo_rem(S_prime, S_minus));
        std::size_t k = 0;
        poly Z;
        std::vector<poly> Ak;
        while((Z = Y - S_star.diff()).is_not_zero()){
            Ak.push_back(gcd(S_star, Z));
            S_star = std::get<0>(quo_rem(S_star, Ak[k]));
            Y = std::get<0>(quo_rem(Z, Ak[k]));
            ++k;
        }
        Ak.push_back(S_star);
        Ak[0] *= c;
        return std::move(Ak);
    }

    static std::vector<poly> squarefree(const poly &A){
        return squarefree_musser(A);
    }

    static aux::vector_map<poly, int> squarefree_factor_list(const poly &A){
        coefficient_type c = A.content();
        poly S = std::get<0>(quo_rem(A, poly(c)));
        poly S_minus = gcd(S, S.diff());
        poly S_star = std::get<0>(quo_rem(S, S_minus));
        int k = 0;
        aux::vector_map<poly, int> inverseAk;
        while(S_minus.deg() > 0){
            poly Y = gcd(S_star, S_minus);
            poly An = std::get<0>(quo_rem(S_star, Y));
            if(An != poly(1)){
                auto p = inverseAk.insert(std::make_pair(An, k + 1));
            }
            S_star = Y;
            S_minus = std::get<0>(quo_rem(S_minus, Y));
            ++k;
        }
        if(S_star != poly(1)){
            auto p = inverseAk.insert(std::make_pair(S_star, k + 1));
        }
        bool find_one = false;
        for(std::pair<poly, int> &i : inverseAk){
            if(i.second == 1){
                i.first *= poly(c) * S_minus;
                find_one = true;
                break;
            }
        }
        if(!find_one){
            inverseAk.insert(std::make_pair(poly(c) * S_minus, 1));
        }
        return std::move(inverseAk);
    }

    template<class T>
    static T abs(const T &a){
        return aux::abs(a);
    }

    template<class T>
    static T pow(T a, degree_type n){
        aux::pow<T>(a, n);
    }

    template<>
    static poly pow<poly>(poly a, degree_type n){
        if(n == 0){
            return poly(1);
        }
        bool positive = n > 0;
        if(!positive){
            n = -n;
        }
        poly r = poly(1);
        for(; n; n >>= 1, a *= a){
            if(n & 1){
                r *= a;
            }
        }
        if(positive){
            return std::move(r);
        }else{
            throw;
        }
    }

    static int gcd(int a, int b){
        a = abs(a);
        b = abs(b);
        if(a < b){
            std::swap(a, b);
        }
        while(b != 0){
            int q = a / b;
            int r = a % b;
            a = b;
            b = r;
        }
        return a;
    }

    static coefficient_type gcd(const coefficient_type &a_, const coefficient_type &b_){
        int a = coe_to_int(a_), b = coe_to_int(b_);
        return gcd(a, b);
    }

    template<class T, class... U>
    static std::vector<T> unpack(const T &a, const U&... seq){
        return std::vector<T>({ a, seq... });
    }

    class coe_to_int_exception : public std::runtime_error{
    public:
        coe_to_int_exception() : std::runtime_error("coefficient's denominator is not 1."){}
        coe_to_int_exception(const coe_to_int_exception&) = default;
    };

    static int coe_to_int(const coefficient_type &coe){
        mpz_class n = coe.get_num();
        mpz_class d = coe.get_den();
        if(d != 1){
            throw coe_to_int_exception();
        }
        return static_cast<int>(n.get_si());
    }

    template<class T>
    static bool odd(const T &a){
        return (a & 1) == 1;
    }

    template<>
    static bool odd<poly>(const poly &a){
        return (coe_to_int(a(0)) & 1) == 1;
    }

    template<class T>
    static bool even(const T &a){
        return (a & 1) == 0;
    }

    template<>
    static bool even<poly>(const poly &a){
        return (coe_to_int(a(0)) & 1) == 0;
    }

public:
    static poly parse(const std::string &str){
        return aux::parse(str);
    }

    class parsing_error : public std::runtime_error{
    public:
        parsing_error() : std::runtime_error("poly: parsing error"){}
        parsing_error(const parsing_error&) = default;
        parsing_error(parsing_error&&) = default;
    };
};

namespace aux{
    enum class token_id : int{
        quo_rem = 3,
        add = 4,
        sub = 5,
        pow = 2,
        x = 6,
        value = 7,
        space = 0,
        end = 2147483647
    };

    class semantic_data{
    public:
        virtual ~semantic_data() = default;
    };

    template<class Iter>
    struct poly_lexer{
        using iterator = Iter;

        struct token_type{
            using identifier_type = token_id;
            token_type(){}
            token_type(const token_type&) = delete;
            token_type(token_type &&other) :
                first(std::move(other.first)),
                last(std::move(other.last)),
                identifier(other.identifier)
            {}

            ~token_type() = default;
            iterator first, last;
            identifier_type identifier;
        };

        template<class Action>
        static std::vector<token_type> tokenize(iterator iter, iterator end, Action &action){
            std::vector<token_type> result;
            iterator first = iter;
            std::size_t line_num = 0, char_num = 0, word_num = 0;
            char c;

            state_1:;
            if(iter == end){
                goto end_of_tokenize;
            }
            c = *iter;
            switch(c){
            case  32: 
                ++char_num;
                ++iter;
                goto state_2;
            case  43: 
                ++char_num;
                ++iter;
                goto state_3;
            case  45: 
                ++char_num;
                ++iter;
                goto state_4;
            case  47: 
                ++char_num;
                ++iter;
                goto state_5;
            case  48: case  49: case  50: case  51: case  52: case  53: case  54: case  55:
            case  56: case  57: 
                ++char_num;
                ++iter;
                goto state_6;
            case  94: 
                ++char_num;
                ++iter;
                goto state_7;
            case 120: 
                ++char_num;
                ++iter;
                goto state_8;
            }
            throw std::runtime_error("lexical error : state 1");

            state_2:;
            if(iter == end){
                goto end_of_tokenize;
            }
            c = *iter;
            switch(c){
            case  32: 
                ++char_num;
                ++iter;
                goto state_2;
            }
            {
                first = iter;
                goto state_1;
            }

            state_3:;
            if(iter == end){
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::add;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::add;
                result.push_back(std::move(t));
                first = iter;
                goto state_1;
            }

            state_4:;
            if(iter == end){
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::sub;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::sub;
                result.push_back(std::move(t));
                first = iter;
                goto state_1;
            }

            state_5:;
            if(iter == end){
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::quo_rem;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::quo_rem;
                result.push_back(std::move(t));
                first = iter;
                goto state_1;
            }

            state_6:;
            if(iter == end){
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::value;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            c = *iter;
            switch(c){
            case  48: case  49: case  50: case  51: case  52: case  53: case  54: case  55:
            case  56: case  57: 
                ++char_num;
                ++iter;
                goto state_6;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::value;
                result.push_back(std::move(t));
                first = iter;
                goto state_1;
            }

            state_7:;
            if(iter == end){
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::pow;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::pow;
                result.push_back(std::move(t));
                first = iter;
                goto state_1;
            }

            state_8:;
            if(iter == end){
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::x;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.identifier = token_type::identifier_type::x;
                result.push_back(std::move(t));
                first = iter;
                goto state_1;
            }

            end_of_tokenize:;
            {
                token_type t;
                t.first = iter;
                t.last = iter;
                t.identifier = token_type::identifier_type::end;
                result.push_back(std::move(t));
            }
            return result;
        }
    };

    struct parsing_semantic_data : public semantic_data{
        parsing_semantic_data() = default;
        parsing_semantic_data(const parsing_semantic_data &other) : p(other.p){}
        parsing_semantic_data(const poly &p) : p(p){}
        poly p;

        template<class Iter>
        semantic_data *make_value(Iter first, Iter last){
            std::string str(first, last);
            Coefficient coe = 0;
            for(auto &i : str){
                coe *= 10;
                coe += i - '0';
            }
            return new parsing_semantic_data(poly(0, coe));
        }
    };

    template<class Lexer>
    class parser{
    public:
        template<class Iter>
        static void parse(poly &poly, Iter iter, Iter end){
            bool first_term = true;
            do{
                iter = term(poly, iter, first_term);
                first_term = false;
            }while(iter->identifier != token_id::end);
        }

        template<class Iter>
        static Iter term(poly &p, Iter iter, bool first_term){
            bool positive = true;
            if(iter->identifier == token_id::add){
                positive = true;
                ++iter;
            }else if(iter->identifier == token_id::sub){
                positive = false;
                ++iter;
            }else if(!first_term){
                throw;
            }
            auto make_str = [](decltype(Iter::value_type::first) beg, decltype(Iter::value_type::last) end){
                std::string str(beg, end);
                return std::move(str);
            };
            std::string num = "1";
            std::string den = "1";
            if(iter->identifier == token_id::value){
                num = make_str(iter->first, iter->last);
                ++iter;
                if(iter->identifier == token_id::quo_rem){
                    ++iter;
                    if(iter->identifier == token_id::value){
                        den = make_str(iter->first, iter->last);
                        ++iter;
                    }else{
                        throw;
                    }
                }
            }
            std::string exponent = "0";
            if(iter->identifier == token_id::x){
                ++iter;
                if(iter->identifier == token_id::pow){
                    ++iter;
                    if(iter->identifier == token_id::value){
                        exponent = make_str(iter->first, iter->last);
                        ++iter;
                    }else{
                        throw;
                    }
                }else{
                    exponent = "1";
                }
            }
            poly q(
                std::atol(exponent.c_str()),
                poly::coefficient_type(std::atol(num.c_str())) / poly::coefficient_type(std::atol(den.c_str())));
            if(positive){
                p += q;
            }else{
                p -= q;
            }
            return iter;
        }
    };

    poly parse(const std::string &str){
        try{
            parsing_semantic_data sa;
            using lexer = poly_lexer<std::string::const_iterator>;
            std::vector<lexer::token_type> result = lexer::tokenize(str.begin(), str.end(), sa);
            poly poly_data;
            parser<lexer>::parse(poly_data, result.begin(), result.end());
            return poly_data;
        }catch(...){
            throw poly::parsing_error();
        }
    }
}

poly operator +(const poly &lhs, const poly &rhs){
    poly r = lhs;
    r += rhs;
    return r;
}

poly operator -(const poly &lhs, const poly &rhs){
    poly r = lhs;
    r -= rhs;
    return r;
}

poly operator *(const poly &lhs, const poly &rhs){
    return poly::mul(lhs, rhs);
}

bool operator <(const poly &a, const poly &b){
    return poly::lexicographical_compare(a, b) < 0;
}

std::ostream &operator <<(std::ostream &os, const poly &p){
    bool first = true;
    for(auto iter = p.rbegin(); iter != p.rend(); ++iter){
        auto &i(*iter);
        if(i.coe > 0){
            if(!first){
                os << " + ";
            }
        }else if(i.coe < 0){
            if(first){
                os << "-";
            }else{
                os << " - ";
            }
        }
        std::string
            num_str = std::to_string(mpz_class(abs(i.coe.get_num())).get_si()),
            den_str = std::to_string(mpz_class(abs(i.coe.get_den())).get_si());
        if(den_str == "1"){
            if(num_str == "1" || num_str == "-1"){
                if(i.deg == 0){
                    os << num_str;
                }
            }else{
                os << num_str;
            }
        }else{
            os << (i.coe > 0 ? i.coe : -i.coe);
        }
        if(i.deg > 0){
            if(poly::abs(i.coe) != 1){
                os << " ";
            }
            os << "x";
            if(i.deg > 1){
                os << "^";
                os << i.deg;
            }
        }
        first = false;
    }
    if(first){
        os << "0";
    }
    return os;
}

class poly_frac;
poly_frac operator +(const poly&, const poly_frac&);
poly_frac operator -(const poly&, const poly_frac&);
poly_frac operator *(const poly&, const poly_frac&);
poly_frac operator /(const poly&, const poly_frac&);
poly_frac operator /(const poly&, const poly&);
std::ostream &operator <<(std::ostream&, const poly_frac&);

class poly_frac{
    friend std::ostream &operator <<(std::ostream&, const poly_frac&);

private:
    using factor_list_type = aux::vector_map<poly, int>;

public:
    class transrate_poly_exception : public std::runtime_error{
    public:
        transrate_poly_exception() : std::runtime_error("denominator != 1"){}
        transrate_poly_exception(const transrate_poly_exception&) = default;
    };

    poly_frac() : num({ { poly(1), 0 } }), den({ { poly(1), 1 } }){}
    poly_frac(const poly_frac &other)
        : num(other.num), den(other.den)
    {}
    
    poly_frac(poly_frac &&other)
        : num(std::move(other.num)), den(std::move(other.den))
    {}

    poly_frac(int n) : num({ { poly(n), 1 } }), den({ { poly(1), 1 } }){}
    poly_frac(const poly &p)
        : num(), den()
    { std::tie(num, den) = reduce(p, poly(1)); }
    
    poly_frac(const poly &p, const poly &q)
        : num(), den()
    { std::tie(num, den) = reduce(p, q); }

    poly_frac operator +(const poly_frac &rhs){
        poly_frac r = *this;
        r += rhs;
        return r;
    }

    poly_frac operator -(const poly_frac &rhs){
        poly_frac r = *this;
        r -= rhs;
        return r;
    }

    poly_frac operator *(const poly_frac &rhs){
        poly_frac r = *this;
        r *= rhs;
        return r;
    }

    poly_frac operator /(const poly_frac &rhs){
        poly_frac r = *this;
        r /= rhs;
        return r;
    }

    poly_frac &operator +=(const poly_frac &rhs){
        std::tie(num, den) = reduce(
            expand(num) * expand(rhs.den) + expand(rhs.num) * expand(den),
            expand(den) * expand(rhs.den)
        );
        return *this;
    }

    poly_frac &operator -=(const poly_frac &rhs){
        std::tie(num, den) = reduce(
            expand(num) * expand(rhs.den) - expand(rhs.num) * expand(den),
            expand(den) * expand(rhs.den)
        );
        return *this;
    }

    poly_frac &operator *=(const poly_frac &rhs){
        coefficient_mul(num, rhs.num.begin(), rhs.num.end());
        coefficient_mul(den, rhs.den.begin(), rhs.den.end());
        std::tie(num, den) = reduce(expand(num), expand(den));
        return *this;
    }

    poly_frac &operator /=(const poly_frac &rhs){
        coefficient_mul(num, rhs.den.begin(), rhs.den.end());
        coefficient_mul(den, rhs.num.begin(), rhs.num.end());
        std::tie(num, den) = reduce(expand(num), expand(den));
        return *this;
    }

    poly_frac &operator +=(const poly &rhs){
        std::tie(num, den) = reduce(expand(num) + rhs * expand(den), expand(den));
        return *this;
    }

    poly_frac &operator -=(const poly &rhs){
        std::tie(num, den) = reduce(expand(num) - rhs * expand(den), expand(den));
        return *this;
    }

    poly_frac &operator *=(const poly &rhs){
        factor_list_type f = poly::squarefree_factor_list(rhs);
        num.insert(f.begin(), f.end());
        std::tie(num, den) = reduce(expand(num), expand(den));
        return *this;
    }

    poly_frac &operator /=(const poly &rhs){
        factor_list_type f = poly::squarefree_factor_list(rhs);
        den.insert(f.begin(), f.end());
        std::tie(num, den) = reduce(expand(num), expand(den));
        return *this;
    }

    poly_frac &operator =(const poly_frac &other){
        num = other.num;
        den = other.den;
        return *this;
    }

    poly_frac &operator =(poly_frac &&other){
        num = std::move(other.num);
        den = std::move(other.den);
        return *this;
    }

    bool is_zero() const{
        return num.begin()->first.is_zero();
    }

    bool is_not_zero() const{
        return !is_zero();
    }

    static std::tuple<factor_list_type, factor_list_type> reduce(const poly &p, const poly &q){
        return reduce(poly::squarefree_factor_list(p), poly::squarefree_factor_list(q));
    }

    static std::tuple<factor_list_type, factor_list_type> reduce(factor_list_type fp, factor_list_type fq){
        do{
            bool continue_flag = false;
            for(auto i = fp.begin(); i != fp.end(); ++i){
                auto j = fq.find(i->first);
                if(j != fq.end()){
                    if(j->second <= i->second){
                        i->second -= j->second;
                        fq.erase(j);
                        if(i->second == 0){
                            fp.erase(i);
                        }
                    }else if(i->second < j->second){
                        j->second -= i->second;
                        fp.erase(i);
                    }
                    continue_flag = true;
                }
                if(continue_flag){
                    break;
                }
            }
            if(continue_flag){
                continue;
            }
        }while(false);
        if(fq.empty()){
            fq = factor_list_type({ { 1, 1 } });
        }
        return std::make_tuple(std::move(fp), std::move(fq));
    }

    static poly expand(const factor_list_type &factor_list){
        poly r = 1;
        for(auto &i : factor_list){
            r *= poly::pow(i.first, i.second);
        }
        return std::move(r);
    }

    static void coefficient_mul(
        factor_list_type &a,
        factor_list_type::const_iterator first, factor_list_type::const_iterator last
    ){
        for(; first != last; ++first){
            auto iter = a.find(first->first);
            if(iter != a.end()){
                iter->second += first->second;
            }else{
                a.insert(iter, *first);
            }
        }
    }

    const poly to_poly() const{
        if(expand(den) != 1){
            throw transrate_poly_exception();
        }
        return expand(num);
    }

    const factor_list_type &get_num() const{
        return num;
    }

    const factor_list_type &get_den() const{
        return den;
    }

    bool operator ==(const poly_frac &other) const{
        return expand(num) == expand(other.num) && expand(den) == expand(other.den);
    }

    bool operator !=(const poly_frac &other) const{
        return !(*this == other);
    }

private:
    factor_list_type num, den;
};

poly_frac operator +(const poly &a, const poly_frac &B){
    poly_frac A = a;
    return A + B;
}

poly_frac operator -(const poly &a, const poly_frac &B){
    poly_frac A = a;
    return A - B;
}

poly_frac operator *(const poly &a, const poly_frac &B){
    poly_frac A = a;
    return A * B;
}

poly_frac operator /(const poly &a, const poly_frac &B){
    poly_frac A = a;
    return A / B;
}

poly_frac operator /(const poly &a, const poly &b){
    poly_frac A = a, B = b;
    return A / B;
}

std::ostream &operator <<(std::ostream &os, const poly_frac &p){
    os << poly_frac::expand(p.get_num());
    const poly_frac::factor_list_type &fl = p.get_den();
    if(!(fl.size() == 1 && fl.begin()->first == 1)){
        os << " / " << poly_frac::expand(p.get_den());
    }
    return os;
}

std::tuple<poly_frac, poly_frac> hermite_reduce_original(const poly &A, const poly &D){
    std::vector<poly> Dn = poly::squarefree(D);
    std::vector<poly> P_An;
    {
        std::vector<poly> Dp = Dn;
        for(int i = 0; i < Dp.size(); ++i){
            Dp[i] = poly::pow(Dp[i], i + 1);
        }
        P_An = poly::partial_fraction(A, Dp);
    }
    poly &P = P_An[0];
    poly_frac g = 0;
    poly_frac h = P + P_An[1] / Dn[0];
    for(int k = 2; k - 1 < Dn.size() && Dn[k - 1].deg() > 0; ++k){
        poly V = Dn[k - 1];
        for(int j = k - 1; j >= 1; --j){
            poly B, C;
            std::tie(B, C) = poly::extended_euclidean(V.diff(), V, -poly::div(P_An[k], poly(j)));
            g = g + B / poly::pow(V, j);
            P_An[k] = -j * C - B.diff();
        }
        h = h + P_An[k] / V;
    }
    return std::make_tuple(g, h);
}

std::tuple<poly_frac, poly_frac> hermite_reduce_quadratic(poly A, poly D){
    poly_frac g = 0;
    std::vector<poly> Dn = poly::squarefree(D);
    for(int i = 2; i - 1 < Dn.size() && Dn[i - 1].deg() > 0; ++i){
        poly V = Dn[i - 1];
        poly U = poly::div(D, poly::pow(V, i));
        for(int j = i - 1; j >= 1; --j){
            poly B, C;
            std::tie(B, C) = poly::extended_euclidean(U * V.diff(), V, -poly::div(A, poly(j)));
            g = g + B / poly::pow(V, j);
            A = -j * C - U * B.diff();
        }
        D = U * V;
    }
    return std::make_tuple(g, A / D);
}

std::tuple<poly_frac, poly_frac> hermite_reduce_liner(poly A, poly D){
    poly_frac g = 0;
    poly Dminus = poly::gcd(D, D.diff());
    poly Dstar = poly::div(D, Dminus);
    while(Dminus.deg() > 0){
        poly Dminus2 = poly::gcd(Dminus, Dminus.diff());
        poly Dminus_star = poly::div(Dminus, Dminus2);
        poly B, C;
        std::tie(B, C) = poly::extended_euclidean(poly::div(-Dstar * Dminus.diff(), Dminus), Dminus_star, A);
        A = C - poly::div(B.diff() * Dstar, Dminus_star);
        g = g + B / Dminus;
        Dminus = Dminus2;
    }
    return std::make_tuple(g, A / Dstar);
}

std::tuple<poly_frac, poly_frac> hermite_reduce(const poly &A, const poly &D){
    return hermite_reduce_liner(A, D);
}

template<class Var>
class multivar_poly{
public:
    using term_key_type = aux::vector_map<Var, int>;
    using term_mapped_type = mpq_class;
    using data_type = aux::vector_map<term_key_type, term_mapped_type>;

    multivar_poly() : data(){}
    multivar_poly(const multivar_poly &other) : data(other.data){}
    multivar_poly(multivar_poly &&other) : data(std::move(other.data)){}

    multivar_poly(const Var &var) : data(){
        term_key_type term;
        term.insert(term_key_type::value_type(var, 1));
        data.insert(data_type::value_type(term, 1));
    }

    multivar_poly(int coe) : data(){
        term_key_type term;
        term.insert(term_key_type::value_type(Var(), 1));
        data.insert(data_type::value_type(term, coe));
    }

    multivar_poly(int coe, const Var &var) : data(){
        term_key_type term;
        term.insert(term_key_type::value_type(var, 1));
        data.insert(data_type::value_type(term, coe));
    }

    multivar_poly(int coe, const Var &var, int n) : data(){
        term_key_type term;
        term.insert(term_key_type::value_type(var, n));
        data.insert(data_type::value_type(term, coe));
    }

    multivar_poly &operator +=(const multivar_poly &other){
        for(const data_type::value_type &other_term : other.data){
            auto iter = data.find(other_term.first);
            if(iter != data.end()){
                iter->second += other_term.second;
                if(iter->second == 0){
                    iter = data.erase(iter);
                }
            }else{
                data.insert(other_term);
            }
        }
        return *this;
    }

    multivar_poly &operator -=(const multivar_poly &other){
        for(const data_type::value_type &other_term : other.data){
            auto iter = data.find(other_term.first);
            if(iter != data.end()){
                iter->second -= other_term.second;
                if(iter->second == 0){
                    iter = data.erase(iter);
                }
            }else{
                data_type::value_type other = other_term;
                other.second = -other.second;
                data.insert(other);
            }
        }
        return *this;
    }

    multivar_poly &operator *=(const multivar_poly &other){
        for(data_type::iterator iter = data.begin(); iter != data.end(); ++iter){
            data_type::value_type &term = *iter;
            for(const data_type::value_type &other_term : other.data){
                for(const term_key_type::value_type &other_var : other_term.first){
                    auto jter = term.first.find(other_var.first);
                    if(jter == term.first.end()){
                        term.first.insert(other_var);
                    }else{
                        jter->second += other_var.second;
                        if(jter->second == 0){
                            jter = term.first.erase(jter);
                        }
                    }
                }
                term.second *= other_term.second;
            }
        }
        return *this;
    }

    multivar_poly &operator /=(const multivar_poly &other){
        for(data_type::iterator iter = data.begin(); iter != data.end(); ++iter){
            data_type::value_type &term = *iter;
            for(const data_type::value_type &other_term : other.data){
                for(const term_key_type::value_type &other_var : other_term.first){
                    auto jter = term.first.find(other_var.first);
                    if(jter == term.first.end()){
                        term_key_type::value_type other = other_var;
                        other.second = -other.second;
                        term.first.insert(other_var);
                    }else{
                        jter->second -= other_var.second;
                    }
                    if(term.second == 0){
                        iter = data.erase(iter);
                    }
                }
                term.second /= other_term.second;
            }
        }
        return *this;
    }

    multivar_poly operator +(const multivar_poly &other){
        multivar_poly r = *this;
        r += other;
        return r;
    }

    multivar_poly operator -(const multivar_poly &other){
        multivar_poly r = *this;
        r -= other;
        return r;
    }

    multivar_poly operator *(const multivar_poly &other){
        multivar_poly r = *this;
        r *= other;
        return r;
    }

    multivar_poly operator /(const multivar_poly &other){
        multivar_poly r = *this;
        r /= other;
        return r;
    }

    class gaussian_elim_no_solution_exception : public std::runtime_error{
    public:
        gaussian_elim_no_solution_exception() : std::runtime_error("gaussian elim: no solution."){}
        gaussian_elim_no_solution_exception(const gaussian_elim_no_solution_exception&) = default;
    };

    static aux::vector_map<term_key_type, mpq_class> gaussian_elim(const std::vector<std::pair<multivar_poly, mpq_class>> &equations){
        aux::vector_map<term_key_type, mpq_class> r;
        std::vector<term_key_type> multivar_set;
        auto multivar_set_insert = [&](const term_key_type &elem){
            auto iter = std::lower_bound(multivar_set.begin(), multivar_set.end(), elem);
            if(iter == multivar_set.end() || *iter != elem){
                multivar_set.insert(iter, elem);
            }
        };
        for(auto &pair : equations){
            for(auto &term : pair.first.data){
                multivar_set_insert(term.first);
            }
        }
        if(multivar_set.size() != equations.size()){
            throw gaussian_elim_no_solution_exception();
        }
        std::size_t N = equations.size();
        std::vector<std::vector<mpq_class>> M(N);
        for(std::size_t i = 0; i < N; ++i){
            M[i].resize(N + 1);
            M[i][N] = equations[i].second;
            auto iter = multivar_set.begin();
            for(std::size_t j = 0; j < N; ++j, ++iter){
                auto find_result = equations[i].first.data.find(*iter);
                if(find_result != equations[i].first.data.end()){
                    M[i][j] = find_result->second;
                }else{
                    M[i][j] = 0;
                }
            }
        }
        for(std::size_t i = 0; i < N; ++i){
            mpq_class pivot = M[i][i];
            for(std::size_t j = 0; j < N + 1; ++j){
                M[i][j] = (1 / pivot) * M[i][j];
            }
            for(std::size_t k = i + 1; k < N; ++k){
                mpq_class mul = M[k][i];
                for(std::size_t n = i; n < N + 1; ++n){
                    M[k][n] = M[k][n] - mul * M[i][n];
                }
            }
        }
        for(int i = static_cast<int>(N - 1); i > 0; --i){
            for(int k = i - 1; k >= 0; --k){
                mpq_class mul = M[k][i];
                for(int n = i; n < N + 1; ++n){
                    M[k][n] = M[k][n] - mul * M[i][n];
                }
            }
        }
        auto iter = multivar_set.begin();
        for(std::size_t i = 0; i < N; ++i, ++iter){
            r.insert(std::make_pair(*iter, M[i][N]));
        }
        return std::move(r);
    }

private:
    data_type data;
};

namespace aux{
    using type_info = int;
    template<class Dummy>
    class type_info_factory_template{
    public:
        template<class T>
        static type_info get(){
            static type_info storage = get2();
            return storage;
        }

    private:
        static type_info get2(){
            static type_info val = 0;
            return val++;
        }
    };

    using type_info_factory = type_info_factory_template<void>;
}

// 初等関数内で発生する全般的な例外クラス．
class elementary_function_exception : public std::runtime_error{
public:
    elementary_function_exception() = delete;
    elementary_function_exception(const elementary_function_exception&) = default;
    elementary_function_exception(const char *message) : std::runtime_error(message){}
};

// 初等関数．
class elementary_function{
public:
    virtual mpf_class approx(std::size_t prec = 64) = 0;
    virtual std::unique_ptr<elementary_function> clone() const = 0;
    virtual bool equal(const elementary_function&) = 0;
    virtual aux::type_info type() const = 0;
    virtual std::string to_string() const = 0;

    // -------- utility --------
    static mpf_class &ln2(){
        static mpf_class storage;
        return storage;
    }

    static mpf_class &sqrt2(){
        static mpf_class storage;
        return storage;
    }

    static mp_bitcnt_t &current_mp_bitcnt_ln2(){
        static mp_bitcnt_t storage;
        return storage;
    }

    static void calc_ln2(mp_bitcnt_t prec){
        if(prec == current_mp_bitcnt_ln2()){
            return;
        }
        ln2().set_prec(prec);
        sqrt2().set_prec(prec);
        current_mp_bitcnt_ln2() = prec;
        ln2() = 0;
        mpf_class u, u_;
        u.set_prec(current_mp_bitcnt_ln2()), u_.set_prec(current_mp_bitcnt_ln2());
        u = u_ = mpf_class(1) / mpf_class(3);
        int i = 1;
        for(mp_bitcnt_t c = 0; c < current_mp_bitcnt_ln2() / 3 + 1; ++c){
            ln2() += 2 * u / i;
            for(int j = 0; j < 2; ++j){
                u *= u_;
            }
            i += 2;
        }
        sqrt2() = 2;
        sqrt2() = sqrt(sqrt2());
    }
};

// 値．
class direct_value : public elementary_function{
private:
    mpf_class data;
    direct_value(const direct_value &other) : data(other.data){}

public:
    direct_value(double a) : data(a){}
    direct_value(const mpf_class &a) : data(a){}
    direct_value(mpf_class &&a) : data(std::move(a)){}

    mpf_class approx(std::size_t prec = 64) override{
        mpf_class r;
        r = data;
        r.set_prec(prec);
        return r;
    }

    std::unique_ptr<elementary_function> clone() const override{
        return std::unique_ptr<elementary_function>(new direct_value(*this));
    }

    bool equal(const elementary_function &other) override{
        return type() == other.type() && data == static_cast<const direct_value&>(other).data;
    }

    aux::type_info type() const override{
        return aux::type_info_factory::get<decltype(*this)>();
    }

    std::string to_string() const override{
        std::stringstream ss;
        ss << data;
        return ss.str();
    }
};

// 除数が0．
class divide_by_zero_exception : public elementary_function_exception{
public:
    divide_by_zero_exception() : elementary_function_exception("divide by zero exception."){}
    divide_by_zero_exception(const divide_by_zero_exception&) = default;
};

// 四則演算の定義．
#define define_operator(name, op) \
    class name : public elementary_function{ \
    private: \
        std::unique_ptr<elementary_function> lhs, rhs; \
        mutable bool approxed = false; \
        mutable mpf_class value_; \
        void check_divide_by_zero(){ \
            if(std::string(#name) == std::string("div") && rhs == 0){ \
                throw divide_by_zero_exception(); \
            } \
        }\
        name( \
            const elementary_function &lhs, \
            const elementary_function &rhs \
        ) : lhs(lhs.clone()), rhs(rhs.clone()){ check_divide_by_zero(); } \
    public: \
        name(double a, double b) : lhs(new direct_value(a)), rhs(new direct_value(b)){ \
            check_divide_by_zero(); \
        } \
        name(const mpf_class &a, const mpf_class &b) : lhs(new direct_value(a)), rhs(new direct_value(b)){ \
            check_divide_by_zero(); \
        } \
        mpf_class approx(std::size_t prec = 64) override{ \
            if(!approxed){ \
                mpf_class r = lhs->approx(prec), s = rhs->approx(prec); \
                value_ = r op s; \
                approxed = true; \
                return value_; \
            }else{ \
                return value_; \
            } \
        } \
        std::unique_ptr<elementary_function> clone() const override{ \
            return std::unique_ptr<elementary_function>(new name(*lhs, *rhs)); \
        } \
        bool equal(const elementary_function &other) override{ \
            return \
                type() == other.type() \
                && lhs->equal(*static_cast<const name&>(other).lhs) \
                && rhs->equal(*static_cast<const name&>(other).rhs); \
        } \
        std::string to_string() const override{ \
            return lhs->to_string() + " " + #op + " " + rhs->to_string(); \
        } \
        aux::type_info type() const override{ \
            return aux::type_info_factory::get<decltype(*this)>(); \
        } \
    }

define_operator(add, +);
define_operator(sub, -);
define_operator(mul, *);
define_operator(div, /);

#undef define_operator

// 対数関数の引数が不正．
class log_exception : public elementary_function_exception{
public:
    log_exception() : elementary_function_exception("log_x y: y is invalid argument."){}
    log_exception(const log_exception&) = default;
};

// 対数関数．
class ln : public elementary_function{
private:
    ln(const ln &other) : argument(nullptr){
        argument = std::move(argument->clone());
    }

public:
    ln() = delete;
    ln(const elementary_function &argument) : argument(argument.clone()){}
    ln(std::unique_ptr<elementary_function> argument) : argument(std::move(argument)){}

    mpf_class approx(std::size_t prec = 64) override{
        if(!approxed){
            mpf_class x = argument->approx(prec);
            if(x <= 0){
                throw log_exception();
            }
            bool inverse = x < 1;
            if(inverse){
                x = 1 / x;
            }
            mp_bitcnt_t prec_bit = x.get_prec();
            calc_ln2(prec_bit);
            mpf_class s;
            s.set_prec(x.get_prec());
            int k = mpf_class(x / sqrt2()).r2_exp();
            {
                mpz_class k_ = 1;
                k_ <<= k;
                x /= k_;
            }
            x -= 1;
            mp_bitcnt_t n = (prec_bit / 32 + 1) * 5;
            for(mp_bitcnt_t i = n; i >= 1; --i){
                s = i * x / (2 + i * x / (2 * i + 1 + s));
            }
            if(inverse){
                value_ = -(ln2() * k + x / (1 + s));
            }else{
                value_ = ln2() * k + x / (1 + s);
            }
            approxed = true;
        }
        return value_;
    }

    std::unique_ptr<elementary_function> clone() const override{
        return std::unique_ptr<elementary_function>(new ln(*this));
    }

    bool equal(const elementary_function &other) override{
        return other.type() == type() && argument->equal(*static_cast<const ln&>(other).argument);
    }

    aux::type_info type() const override{
        return aux::type_info_factory::get<decltype(*this)>();
    }

    std::string to_string() const override{
        return "ln(" + argument->to_string() + ")";
    }

private:
    std::unique_ptr<elementary_function> argument;
    mutable bool approxed = false;
    mutable mpf_class value_;
};

// 指数関数．
class exp : public elementary_function{
private:
    exp(const exp &other) : argument(nullptr){
        argument = std::move(argument->clone());
    }

public:
    exp() = delete;
    exp(const elementary_function &argument) : argument(argument.clone()){}
    exp(std::unique_ptr<elementary_function> argument) : argument(std::move(argument)){}

    mpf_class approx(std::size_t prec = 64) override{
        if(!approxed){
            mpf_class x = argument->approx(prec);
            if(x == 0){ return mpf_class(1); }
            bool inverse = x < 0;
            if(inverse){
                x = -x;
            }
            mp_bitcnt_t prec_bit = x.get_prec(), n = (prec_bit / 128 + 1) * 50;
            calc_ln2(prec_bit);
            mpz_class k;
            mpf_class x2, w, tmp, half;
            x2.set_prec(prec_bit), w.set_prec(prec_bit), tmp.set_prec(prec_bit), half.set_prec(prec_bit);
            k = x / ln2() + (x >= 0 ? 0.5 : -0.5);
            x -= k * ln2();
            x2 = x * x;
            w = x2 / n;
            for(mp_bitcnt_t i = n - 4; i >= 6; i -= 4){
                w = x2 / (w + i);
            }
            tmp = (2 + w + x) / (2 + w - x);
            mp_bitcnt_t
                l = mpz_class(k / (sizeof(mp_bitcnt_t) * 8)).get_ui(),
                m = mpz_class(k % (sizeof(mp_bitcnt_t) * 8)).get_ui();
            if(l != 0 || m != 0){
                tmp.dexp(static_cast<mp_exp_t>(l), static_cast<mp_exp_t>(m));
                value_ = inverse ? 1 / tmp : tmp;
            }else{
                value_ = 1 / tmp;
            }
            approxed = true;
        }
        return value_;
    }

    virtual std::unique_ptr<elementary_function> clone() const override{
        return std::unique_ptr<elementary_function>(new exp(*this));
    }

    bool equal(const elementary_function &other) override{
        return other.type() == type() && argument->equal(*static_cast<const exp&>(other).argument);
    }

    aux::type_info type() const override{
        return aux::type_info_factory::get<decltype(*this)>();
    }

    std::string to_string() const override{
        return "exp(" + argument->to_string() + ")";
    }

private:
    std::unique_ptr<elementary_function> argument;
    mutable bool approxed = false;
    mutable mpf_class value_;
};

}

//-------- include end --------//

template<class Container>
Container standard_suffle(Container c){
    Container r;
    std::mt19937 mt(0xFF78FF);
    std::size_t n = c.size();
    for(std::size_t i = 0; i < n; ++i){
        std::size_t m = n - i;
        std::size_t k = mt() % m;
        r.insert(r.end(), std::move(c[k]));
        c.erase(c.begin() + k);
    }
    return std::move(r);
}

#include <cmath>

void elementary_function_test(){
    symbolic_alg::direct_value v(50);
    symbolic_alg::exp e(v);
    std::cout << e.approx(64).get_str() << std::endl;

    symbolic_alg::ln ln(v);
    std::cout << ln.approx(64).get_str() << std::endl;
}

void squarefree_test(){
    using namespace symbolic_alg;
    auto result = poly::squarefree_factor_list(poly::parse("x^5 - x^3 - x^2 + 1"));
    for(auto &i : result){
        std::cout << i.first << " : " << i.second << std::endl;
    }
}

void poly_frac_test(){
    using namespace symbolic_alg;
    poly_frac pf = poly::parse("x^5 - x^3 - x^2 + 1") / poly::parse("x^3 + 2x^2 + 2x + 1") / poly::parse("-x + 1");
    std::cout << pf << std::endl;
    pf *= poly::parse("x^3 + 2x^2 + 2x + 1") * poly::parse("-x + 1");
    std::cout << pf << std::endl;
}

void hermite_reduce_test(){
    using namespace symbolic_alg;
    {
        std::cout << "Hermite Reduce (original version)." << std::endl;
        auto p = hermite_reduce_original(poly::parse("x^7 - 24x^4 - 4x^2 + 8x - 8"), poly::parse("x^8 + 6x^6 + 12x^4 + 8x^2"));
        std::cout << std::get<0>(p) << std::endl;
        std::cout << std::get<1>(p) << std::endl << std::endl;
    }

    {
        std::cout << "Hermite Reduce (quadratic version)." << std::endl;
        auto p = hermite_reduce_quadratic(poly::parse("x^7 - 24x^4 - 4x^2 + 8x - 8"), poly::parse("x^8 + 6x^6 + 12x^4 + 8x^2"));
        std::cout << std::get<0>(p) << std::endl;
        std::cout << std::get<1>(p) << std::endl << std::endl;
    }

    {
        std::cout << "Hermite Reduce (liner version)." << std::endl;
        auto p = hermite_reduce_liner(poly::parse("x^7 - 24x^4 - 4x^2 + 8x - 8"), poly::parse("x^8 + 6x^6 + 12x^4 + 8x^2"));
        std::cout << std::get<0>(p) << std::endl;
        std::cout << std::get<1>(p) << std::endl << std::endl;
    }
}

void multivar_poly_test_1(){
    using namespace symbolic_alg;
    using mvpoly = multivar_poly<std::string>;
    mvpoly a(2, "x", 2);
    mvpoly b(3, "x");
    mvpoly c(4, "x", 3);
    mvpoly d(1, "y", 3);
    mvpoly e(1, "x", -6);
    mvpoly f(1, "y", -3);
    a *= b;
    a *= c;
    a *= d;
    a *= e;
    a *= f;
    return;
}

void multivar_poly_test_2(){
    using namespace symbolic_alg;
    using mvpoly = multivar_poly<std::string>;
    mvpoly a(2, "x", 2);
    mvpoly b(3, "x");
    mvpoly c(1, "x", 2);
    mvpoly d(3, "x");
    mvpoly e(3, "x", 2);
    a += b;
    a += c;
    a += d;
    a -= e;
    return;
}

void gaussian_elim_test(){
    using namespace symbolic_alg;
    using mvpoly = multivar_poly<std::string>;
    std::vector<std::pair<mvpoly, mpq_class>> equations;
    {
        mvpoly lhs;
        lhs += mvpoly(1, "x");
        lhs += mvpoly(1, "y");
        lhs += mvpoly(1, "z");
        lhs += mvpoly(1, "w");
        mpq_class rhs = 10;
        equations.push_back(std::make_pair(lhs, rhs));
    }
    {
        mvpoly lhs;
        lhs += mvpoly(2, "x");
        lhs += mvpoly(1, "y");
        lhs += mvpoly(2, "z");
        lhs += mvpoly(1, "w");
        mpq_class rhs = 14;
        equations.push_back(std::make_pair(lhs, rhs));
    }
    {
        mvpoly lhs;
        lhs += mvpoly(1, "x");
        lhs += mvpoly(2, "y");
        lhs += mvpoly(3, "z");
        lhs += mvpoly(-4, "w");
        mpq_class rhs = -2;
        equations.push_back(std::make_pair(lhs, rhs));
    }
    {
        mvpoly lhs;
        lhs += mvpoly(1, "x");
        lhs += mvpoly(-1, "y");
        lhs += mvpoly(-1, "z");
        lhs += mvpoly(1, "w");
        mpq_class rhs = 0;
        equations.push_back(std::make_pair(lhs, rhs));
    }
    auto solutions = mvpoly::gaussian_elim(equations);
    return;
}

int main(){
    gaussian_elim_test();
    return 0;
}
