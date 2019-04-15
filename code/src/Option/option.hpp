#ifndef OPTION_H
#define OPTION_H

#include <functional>

class Option{

    public:

        Option(double _horizon, double _strike, std::function<double (double) > _payoff);
        ~Option(){};
        std::function<double (double) > payoff;
        // Fonctions inutiles pour Monte Carlo, uniquement pour EDP
        virtual double boundary_left(double t, double x, double rate){ return 0;};
        virtual double boundary_right(double t, double x, double rate){ return 0;};

        double Horizon() const { return horizon;};
        double Strike() const { return strike;}
    private:
        double horizon;
        double strike;
};

#endif