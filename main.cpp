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
#include <boost/rational.hpp>

template<class DegreeType, class CoeffidientType>
class poly;

namespace aux_poly{
    template<class Degree, class Coefficient>
    struct parsing_semantic_data;

    template<class Degree, class Coefficient>
    poly<Degree, Coefficient> parse(const std::string &str);
}

template<class DegreeType, class CoeffidientType>
class poly{
    friend aux_poly::parsing_semantic_data<DegreeType, CoeffidientType>;

public:
    // éüêîÇÃå^
    using degree_type = DegreeType;

    // åWêîÇÃå^
    using coefficient_type = CoeffidientType;

    // çÄÇÃå^
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

    static std::tuple<poly, poly> euclidean_div_int(const poly &A, const poly &B){
        poly Q = 0, R = A;
        degree_type delta;
        while(R.is_not_zero() && (delta = R.deg() - B.deg()) >= 0){
            poly T(delta, coe_to_int(R.lc()) / coe_to_int(B.lc()));
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

    static std::tuple<poly, poly> div(const poly &A, const poly &B){
        return euclidean_div(A, B);
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
        for(const term_type &i : data){
            for(const term_type &j : other.data){
                if(i.deg != j.deg || i.coe != j.coe){
                    return false;
                }
            }
        }
        return true;
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
            if(t.coe.denominator() != 1){
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
            int n = std::abs(data[0].coe.numerator()), d = std::abs(data[0].coe.denominator());
            for(std::size_t i = 1; i < data.size(); ++i){
                if(!(data[i].coe.numerator() % n == 0 && data[i].coe.denominator() % d == 0)){
                    bool divid = false;
                    if(n % data[i].coe.numerator() == 0){
                        n = std::abs(data[i].coe.numerator());
                        divid = true;
                    }
                    if(d % data[i].coe.denominator() == 0){
                        d = std::abs(data[i].coe.denominator());
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
            auto qr = div(a, b);
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
        auto tr = div(std::get<1>(sg) - std::get<0>(sg) * a, b);
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
        std::tie(q, r) = div(c, g);
        if(r.is_not_zero()){
            throw extended_euclidean_exception("c is not in the ideal generated by a and b.");
        }
        s = q * s;
        t = q * t;
        if(s.is_not_zero() && s.deg() >= b.deg()){
            std::tie(q, r) = div(s, b);
            s = r;
            t = t + q * a;
        }
        return std::make_tuple(s, t);
    }

    static poly half_extended_euclidean(const poly &a, const poly &b, const poly &c){
        poly s, g, q, r;
        std::tie(s, g) = half_extended_euclidean(a, b);
        std::tie(q, r) = div(c, g);
        if(r.is_not_zero()){
            throw extended_euclidean_exception("c is not in the ideal generated by a and b.");
        }
        s = q * s;
        if(s.is_not_zero() && s.deg() >= b.deg()){
            std::tie(q, r) = div(s, b);
            s = r;
        }
        return s;
    }

    static std::tuple<poly, poly> half_full_extended_euclidean(const poly &a, const poly &b, const poly &c){
        poly s = half_extended_euclidean(a, b, c), t, r;
        std::tie(t, r) = div(c - s * a, b);
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
        std::tie(a0, r) = div(a, d_prod);
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
                std::tie(q, aij) = div(an[i], d[i]);
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
        std::vector<coefficient_type> É¡, É¿, r;
        std::vector<degree_type> É¬;
        É¡.push_back(0); É¡.push_back(-1);
        É¬.push_back(0); É¬.push_back(A.deg() - B.deg());
        É¿.push_back(0); É¿.push_back(pow(-1, É¬[1] + 1));
        r.push_back(0);
        while(R[i].is_not_zero()){
            r.push_back(R[i].lc());
            poly Quo, Rem;
            std::tie(Quo, Rem) = psude_div(R[i - 1], R[i]);
            R.push_back(std::get<0>(div(Rem, poly(É¿[i]))));
            ++i;
            É¡.push_back(pow(-É¡[i - 1], É¬[i - 1]) * pow(É¡[i - 1], 1 - É¬[i - 1]));
            É¬.push_back(R[i - 1].deg() - R[i].deg());
            É¿.push_back(-r[i - 1] * pow(É¡[i], É¬[i]));
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
            c = c * pow(std::get<0>(div(É¿[j], pow(r[j], 1 + É¬[j]))), R[j].deg()) * pow(r[j], R[j - 1].deg() - R[j + 1].deg());
        }
        return std::make_tuple(poly(s) * poly(c) * pow(R[k], R[k - 1].deg()), std::move(R));
    }

    static std::vector<poly> squarefree_musser(const poly &A){
        coefficient_type c = A.content();
        poly S = std::get<0>(div(A, poly(c)));
        poly S_minus = gcd(S, S.diff());
        poly S_star = std::get<0>(div(S, S_minus));
        std::size_t k = 0;
        std::vector<poly> Ak;
        while(S_minus.deg() > 0){
            poly Y = gcd(S_star, S_minus);
            Ak.push_back(std::get<0>(div(S_star, Y)));
            S_star = Y;
            S_minus = std::get<0>(div(S_minus, Y));
            ++k;
        }
        Ak.push_back(S_star);
        Ak[0] *= poly(c) * S_minus;
        return std::move(Ak);
    }

    static std::vector<poly> squarefree_yun(const poly &A){
        coefficient_type c = A.content();
        poly S = std::get<0>(div(A, poly(c)));
        poly S_prime = S.diff();
        poly S_minus = gcd(S, S_prime);
        poly S_star = std::get<0>(div(S, S_minus));
        poly Y = std::get<0>(div(S_prime, S_minus));
        std::size_t k = 0;
        poly Z;
        std::vector<poly> Ak;
        while((Z = Y - S_star.diff()).is_not_zero()){
            Ak.push_back(gcd(S_star, Z));
            S_star = std::get<0>(div(S_star, Ak[k]));
            Y = std::get<0>(div(Z, Ak[k]));
            ++k;
        }
        Ak.push_back(S_star);
        Ak[0] *= c;
        return std::move(Ak);
    }

    static std::vector<poly> squarefree(const poly &A){
        return squarefree_musser(A);
    }

    //static std::tuple<poly, poly> hermite_reduce(poly A, const poly &D){
    //    poly g;
    //    poly D_minus = gcd(D, D.diff());
    //    poly D_star = D / D_minus;
    //    while(D_minus.deg() > 0){
    //        poly D_minus2 = gcd(D_minus, D_minus.diff());
    //        poly D_minus_star = D_minus / D_minus2;
    //        poly B, C;
    //        std::tie(B, C) = extended_euclidean(-D_star * D_minus.diff() / D_minus, D_minus_star, A);
    //        A = C - B.diff() * D_star / D_minus_star;
    //        g = g + B / D_minus;
    //        D_minus = D_minus2;
    //    }
    //    return std::make_tuple(g, std::get<0>(div(A, D_star)));
    //}

    template<class T>
    static T abs(const T &a){
        return a < 0 ? -a : a;
    }

    template<class T>
    static T pow(T a, degree_type n){
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
        a = std::abs(a);
        b = std::abs(b);
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
        int n = coe.numerator();
        int d = coe.denominator();
        if(d != 1){
            throw coe_to_int_exception();
        }
        return n;
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
        return aux_poly::parse<DegreeType, CoeffidientType>(str);
    }

    class parsing_error : public std::runtime_error{
    public:
        parsing_error() : std::runtime_error("poly: parsing error"){}
        parsing_error(const parsing_error&) = default;
        parsing_error(parsing_error&&) = default;
    };
};

namespace aux_poly{
    enum class token_id : int{
        div = 3,
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
                line_num(other.line_num),
                char_num(other.char_num),
                word_num(other.word_num),
                identifier(other.identifier)
            {}

            ~token_type() = default;
            iterator first, last;
            std::size_t line_num, char_num, word_num;
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
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
                t.identifier = token_type::identifier_type::add;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
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
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
                t.identifier = token_type::identifier_type::sub;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
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
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
                t.identifier = token_type::identifier_type::div;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
                t.identifier = token_type::identifier_type::div;
                result.push_back(std::move(t));
                first = iter;
                goto state_1;
            }

            state_6:;
            if(iter == end){
                token_type t;
                t.first = first;
                t.last = iter;
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
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
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
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
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
                t.identifier = token_type::identifier_type::pow;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
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
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
                t.identifier = token_type::identifier_type::x;
                result.push_back(std::move(t));
                goto end_of_tokenize;
            }
            {
                token_type t;
                t.first = first;
                t.last = iter;
                t.line_num = line_num;
                t.char_num = char_num;
                t.word_num = word_num++;
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
                t.line_num = 0;
                t.char_num = 0;
                t.word_num = 0;
                t.identifier = token_type::identifier_type::end;
                result.push_back(std::move(t));
            }
            return result;
        }
    };

    template<class Degree, class Coefficient>
    struct parsing_semantic_data : public semantic_data{
        parsing_semantic_data() = default;
        parsing_semantic_data(const parsing_semantic_data &other) : p(other.p){}
        parsing_semantic_data(const poly<Degree, Coefficient> &p) : p(p){}
        poly<Degree, Coefficient> p;

        template<class Iter>
        semantic_data *make_value(Iter first, Iter last){
            std::string str(first, last);
            Coefficient coe = 0;
            for(auto &i : str){
                coe *= 10;
                coe += i - '0';
            }
            return new parsing_semantic_data(poly<Degree, Coefficient>(0, coe));
        }
    };

    template<class Lexer>
    class parser{
    public:
        template<class Degree, class Coefficient, class Iter>
        static void parse(poly<Degree, Coefficient> &poly, Iter iter, Iter end){
            bool first_term = true;
            do{
                iter = term(poly, iter, first_term);
                first_term = false;
            }while(iter->identifier != token_id::end);
        }

        template<class Degree, class Coefficient, class Iter>
        static Iter term(poly<Degree, Coefficient> &p, Iter iter, bool first_term){
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
                if(iter->identifier == token_id::div){
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
            poly<Degree, Coefficient> q(
                Degree(std::atol(exponent.c_str())),
                Coefficient(std::atol(num.c_str())) / Coefficient(std::atol(den.c_str()))
            );
            if(positive){
                p += q;
            }else{
                p -= q;
            }
            return iter;
        }
    };

    template<class Degree, class Coefficient>
    poly<Degree, Coefficient> parse(const std::string &str){
        try{
            parsing_semantic_data<Degree, Coefficient> sa;
            using lexer = poly_lexer<std::string::const_iterator>;
            std::vector<lexer::token_type> result = lexer::tokenize(str.begin(), str.end(), sa);
            poly<Degree, Coefficient> poly_data;
            parser<lexer>::parse(poly_data, result.begin(), result.end());
            return poly_data;
        }catch(...){
            throw poly<Degree, Coefficient>::parsing_error();
        }
    }
}

template<class Degree, class Coefficient>
poly<Degree, Coefficient> operator +(const poly<Degree, Coefficient> &lhs, const poly<Degree, Coefficient> &rhs){
    poly<Degree, Coefficient> r = lhs;
    r += rhs;
    return r;
}

template<class Degree, class Coefficient>
poly<Degree, Coefficient> operator -(const poly<Degree, Coefficient> &lhs, const poly<Degree, Coefficient> &rhs){
    poly<Degree, Coefficient> r = lhs;
    r -= rhs;
    return r;
}

template<class Degree, class Coefficient>
poly<Degree, Coefficient> operator *(const poly<Degree, Coefficient> &lhs, const poly<Degree, Coefficient> &rhs){
    return poly<Degree, Coefficient>::mul(lhs, rhs);
}

template<class Degree, class Coefficient>
std::ostream &operator <<(std::ostream &os, const poly<Degree, Coefficient> &p){
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
        std::stringstream ss;
        ss << (i.coe > 0 ? i.coe : -i.coe);
        std::string coe_str = ss.str();
        std::string num_str, den_str;
        auto jter = coe_str.begin();
        while(jter != coe_str.end() && *jter != '/'){
            num_str += *jter;
            ++jter;
        }
        if(num_str == "0"){
            den_str = "1";
        }else{
            if(jter != coe_str.end()){
                ++jter;
                while(jter != coe_str.end()){
                    den_str += *jter;
                    ++jter;
                }
            }
        }
        if(den_str == "1"){
            if(num_str == "1"){
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
            if(poly<Degree, Coefficient>::abs(i.coe) != 1){
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

int main(){
    using rational = boost::rational<int>;
    using p = poly<int, rational>;
    //{
    //    std::cout << "Hermite Reduce\n";
    //    auto r = p::hermite_reduce(p::parse("x^7 - 24x^4 - 4x^2 + 8x - 8"), p::parse("x^8 + 6x^6 + 12x^4 + 8x^2"));
    //    std::cout << std::get<0>(r) << std::endl;
    //    std::cout << std::get<1>(r) << std::endl;
    //}

    {
        std::cout << "Squarefree Musser's algorithm\n";
        auto r = p::squarefree_musser(p::parse("x^8 + 6x^6 + 12x^4 + 8x^2"));
        for(auto &&i : r){
            std::cout << i << "\n";
        }
    }

    {
        std::cout << "Squarefree Yun's algorithm\n";
        auto r = p::squarefree_yun(p::parse("x^8 + 6x^6 + 12x^4 + 8x^2"));
        for(auto &&i : r){
            std::cout << i << "\n";
        }
    }

    {
        std::cout << "Squarefree test" << std::endl;
        auto r = p::squarefree(p::parse("x^5 - x^3 - x^2 + 1"));
        for(auto &&i : r){
            std::cout << i << "\n";
        }
    }

    return 0;
}
