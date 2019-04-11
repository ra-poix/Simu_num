#include <iostream>
#include <cmath>
#include "compose.hpp"
#include "Call.hpp"
#include "BSmodel.hpp"
#include <random>


struct mean_var {
    mean_var(unsigned n = 0, double sum_x = 0, double sum_xx = 0)
    : sample_size(n), sum_x(sum_x), sum_xx(sum_xx) { }
    double mean() const { return sum_x / (double) sample_size; }
    double var() const { return (sum_xx - sample_size * mean() * mean())
        / (double) (sample_size-1); }
    double ic_size() const { return 1.96 * std::sqrt(var() / sample_size); }
    
    mean_var & operator+=(mean_var const & mv) {
        sample_size += mv.sample_size;
        sum_x += mv.sum_x;
        sum_xx += mv.sum_xx;
        return *this;
    }
    friend mean_var operator+(mean_var const & mv1, mean_var const & mv2) {
        return { mv1.sample_size + mv2.sample_size,
            mv1.sum_x + mv2.sum_x,
            mv1.sum_xx + mv2.sum_xx };
    }
    friend mean_var operator*(double alpha, mean_var const & mv) {
        return { mv.sample_size, alpha * mv.sum_x, alpha * mv.sum_xx };
    }
    friend std::ostream & operator<<(std::ostream & o, mean_var const & mv) {
        return o << "Size: " << mv.sample_size
        << "\tMean: " << mv.mean()
        << "\tVar: " << mv.var();
    }
protected:
    unsigned sample_size;
    double sum_x, sum_xx;
};

template <typename TDistrib, typename TGen>
mean_var monte_carlo(TDistrib & X, TGen & gen, unsigned sample_size) {
    double sum_x = 0, sum_xx = 0;
    for (unsigned k = 0; k < sample_size; ++k) {
        double x = X(gen);
        sum_x += x;
        sum_xx += x*x;
    }
    return { sample_size, sum_x, sum_xx };
}

template <typename TDistrib, typename TGen>
mean_var monte_carlo(TDistrib & X, TGen & gen, unsigned batch_size, double epsilon) {
    auto r = monte_carlo(X, gen, batch_size);
    while (r.ic_size() > epsilon) {
        auto tmp = monte_carlo(X, gen, batch_size);
        r += tmp;
    }
    return r;
}


struct sensib_malliavin{
public:
    sensib_malliavin(double delta, double gamma, double vega): delta(0), gamma(0), vega(0){}
    ~sensib_malliavin(){};
protected:
    double delta;
    double gamma;
    double vega;
};

template <typename TDistrib, typename TGen>
sensib_malliavin malliavin(TDistrib & X, TGen & gen, unsigned batch_size, double epsilon){
    
}

template <typename TDistrib, typename TGen>
sensib_malliavin malliavin(composed<Call ,BSmodel>, TGen & gen, unsigned batch_size, double epsilon){
    std::normal_distribution<> G;
    auto Y = compose(call, G);
    auto price = monte_carlo(Y, gen, 1e3, 1e-3);
}

template <typename TDistrib, typename TGen>
auto delta_call_bs(Call x, BSmodel y, TDistrib G ){
    double S = exp((y.Rate(0.0,0.0)-0.5*pow(y.Sigma(0.0,0.0),2))*x.Horizon()+G);
    return exp(-y.Rate(0.0,0.0))*x.Horizon())*call::payoff(exp((y.Rate(0.0,0.0)+));
};
